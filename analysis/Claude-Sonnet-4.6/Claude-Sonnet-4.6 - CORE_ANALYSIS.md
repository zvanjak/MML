# MML Core Layer & Interfaces Analysis

**Analyzed by:** Claude Sonnet 4.6  
**Date:** 2026-03-13  
**Scope:** `mml/interfaces/` (all interface files) and `mml/core/` (all subdirectories)

---

## 1. Scope Overview

The core layer is the mathematical engine connecting the base containers to the higher-level algorithms. It implements the most fundamental numerical methods that everything else depends on:

| File / Folder | Purpose |
|---|---|
| `interfaces/IFunction.h` | Complete function interface hierarchy |
| `interfaces/IODESystem*.h` | ODE/DAE system and stepper interfaces |
| `interfaces/IDynamicalSystem.h` | Dynamical system interface |
| `interfaces/ICoordTransf.h` | Coordinate transformation interface |
| `interfaces/ITensor.h`, `ITensorField.h` | Tensor and tensor field interfaces |
| `interfaces/IInterval.h`, `IParametrized.h` | Supporting interfaces |
| `core/Derivation/` | Numerical differentiation (orders 1–8) |
| `core/Integration/` | 1-D to 3-D integration, Gauss quadrature, path, surface, Monte Carlo |
| `core/LinAlgEqSolvers/` | Direct (LU, QR, SVD, Cholesky) and iterative linear solvers |
| `core/CoordTransf/` | Cartesian, spherical, cylindrical, cylindrical, Lorentz transforms |
| `core/OrthogonalBasis/` | Chebyshev, Hermite, Laguerre, Legendre expansion bases |
| `core/FieldOperations.h` | Gradient, divergence, curl, Laplacian in multiple coord. systems |
| `core/MetricTensor.h` | Metric tensor field for general curvilinear coordinates |
| `core/AlgorithmTypes.h` | `AlgorithmStatus` enum and `ConfigBase` pattern |
| `core/RichardsonExtrapolation.h` | Richardson / Aitken extrapolation utilities |
| `core/SingularityHandling.h` | Regularization helpers for near-singular integrals |
| `core/NumericValidation.h` | Input sanity checking utilities |
| `core/FunctionHelpers.h` | Adapters for function types (composition, binding, sampling) |
| `core/MatrixUtils.h` | Norm computation, condition number, matrix analysis |
| `core/CoordSystem.h` | Coordinate system registry |
| `core/Curves.h`, `Surfaces.h` | Named parametric curves and surfaces catalog |
| `core/Jacobians.h` | Jacobian matrix computation |

---

## 2. Strengths

### 2.1 Interface Architecture (mml/interfaces/)

**Rich, deep interface hierarchy:**  
`IFunction.h` defines seven distinct function categories, each a separate interface:

```
IFunction<RetType, ArgType>
    ├── IRealFunction              — f: ℝ → ℝ
    ├── IRealFunctionParametrized  — f: ℝ → ℝ with named parameters
    ├── IScalarFunction<N>         — f: ℝᴺ → ℝ
    ├── IVectorFunction<N>         — f: ℝᴺ → ℝᴺ
    ├── IVectorFunctionNM<N,M>     — f: ℝᴺ → ℝᴹ
    ├── IParametricCurve<N>        — γ: ℝ → ℝᴺ
    └── IParametricSurface<N>      — σ: ℝ² → ℝᴺ
```

This hierarchy covers the complete space of mathematical functions relevant to numerical computing. Every core algorithm (`Derivation`, `Integration`, `FieldOperations`, `PathIntegration`, `SurfaceIntegration`) is templated over these interfaces — creating total separation between algorithm implementation and function representation.

**`IRealFunction::GetValues()` sampling utility:**  
The method samples the function across an equispaced grid and returns `(x, f(x))` pairs in one call. This is a small but practical convenience that avoids boilerplate in visualization and interpolation code.

**ODE interface decomposition:**  
The ODE interface is split across four files:
- `IODESystem.h` — system definition
- `IODESystemDAE.h` — DAE variant
- `IODESystemStepCalculator.h` — single-step compute interface
- `IODESystemStepper.h` — adaptive stepping interface

This separation of concerns allows fixed-step and adaptive integrators to share the same system interface while composing different stepper strategies — a textbook Strategy pattern.

**`IParametrized.h` — generic parameterized object:**  
The `IParametrized` base allows dynamical systems, function families, and physical models to uniformly declare, query, and set named parameters with defined ranges — enabling generic GUI tools, optimization over parameters, and serialization without knowing the specific system type.

---

### 2.2 Derivation Layer (core/Derivation/)

**Analytically optimal step sizes:**  
`DerivationBase.h` derives the step sizes from error analysis:

| Method | Step Formula | Order |
|---|---|---|
| `NDer1` | `h = 2√ε_m` | O(h) |
| `NDer2` | `h = (3ε_m)^(1/3)` | O(h²) |
| `NDer4` | `h = (11.25ε_m)^(1/5)` | O(h⁴) |
| `NDer6` | `h = (...)^(1/7)` | O(h⁶) |
| `NDer8` | `h = (...)^(1/9)` | O(h⁸) |

These step sizes minimize total error (roundoff + truncation) for the respective stencils, derived from classical numerical analysis. This is not guesswork; it is the theoretically correct approach documented with the derivation. Very few libraries make these choices explicit and justified.

**`ScaleStep(h, x)` for large magnitudes:**  
When `|x|` is large, a fixed absolute step `h` is tiny relative to `x` and causes catastrophic cancellation. The function adaptively scales `h` by `|x|` in this regime. This is a common source of bugs in numerical differentiation code; handling it centrally avoids it everywhere.

**Complete derivative coverage across all function types:**  
- `DerivationRealFunction` — f'(x), f''(x), f'''(x)
- `DerivationScalarFunction<N>` — partial derivatives, gradient, Hessian
- `DerivationVectorFunction<N>` — component-wise derivatives
- `DerivationParametricCurve<N>` — curve tangent, normal, binormal (Frenet frame)
- `DerivationParametricSurface<N>` — partial tangent vectors, normal, shape operator
- `DerivationTensorField` — tensor field derivatives
- `Jacobians.h` — numerical Jacobian for vector functions

This completeness means users do not need to implement specialized differentiation for curves or surfaces — the library does it correctly.

**Error estimation:**  
```cpp
Real df = Derivation::NDer4(f, x0, &error);
```
The optional error output parameter is elegant — it does not force every caller to deal with error estimates, but provides them when needed.

---

### 2.3 Integration Layer (core/Integration/)

**`IntegrationResult` with convergence reporting:**  
All integration methods return `IntegrationResult{value, error_estimate, iterations, function_evaluations, converged}`. Users always know whether the integral converged and can programmatically handle failure — a major improvement over bare `double` returns.

**`IntegrateRomberg` — Richardson extrapolation to high order:**  
Romberg builds a triangular Richardson table, doubling the trapezoid count at each step and extrapolating. This achieves exponential convergence for smooth functions with very few function evaluations. The implementation is correct including the convergence check on the full extrapolation table, not just the last two values.

**`TrapIntegrator` — progressive refinement:**  
The trapezoidal integrator tracks the sum of existing function values, adding only new points at each refinement level. This avoids redundant function evaluations across levels — important when the function is expensive.

**`IntegrateGaussKronrod` — adaptive high-accuracy quadrature:**  
The Gauss-Kronrod rules (e.g., G10K21 = 10-point Gauss + 21-point Kronrod) reuse the Gauss points in the Kronrod rule, giving an embedded error estimate with near-no extra cost. The 21-point rule is the recommended default for most production integrals.

**Comprehensive coverage:**  
| Category | Methods Available |
|---|---|
| Basic 1D | Trap, Simpson, Romberg |
| Gaussian | GaussLegendre (multiple orders), GaussHermite, GaussLaguerre, GaussKronrod |
| Improper | IntegrateUpperInf, IntegrateLowerInf, IntegrateInfInf |
| 2D | Standard + adaptive (variable limits) |
| 3D | Standard + adaptive (variable limits via lambdas) |
| Monte Carlo | MonteCarloIntegration with configurable sample count |
| Path | Line integrals (scalar and vector field) |
| Surface | Flux integrals over parametric and body surfaces |

This coverage is on par with scientific computing libraries an order of magnitude larger.

**Variable-limit 2D/3D integration:**  
The API for triangular domains, tetrahedra, and general regions via:
```cpp
auto y_lo = [](Real x) { return 0.0; };
auto y_hi = [](Real x) { return 1.0 - x; };
```
Lambda-based bounds are clean, flexible, and do not require special types. This design handles virtually any integration domain.

---

### 2.4 Linear Algebra Solvers (core/LinAlgEqSolvers/)

**`GaussJordanSolver` — robust pivoting strategy:**  
Full pivoting (both row and column) maximizes numerical stability beyond partial pivoting alone. The singularity threshold uses `ε × ‖A‖∞ × n` (infinity-norm scaled by dimension) — the theoretically correct choice, not a fixed threshold that fails for large or small matrices.

Non-finite pivot detection:
```cpp
if constexpr (std::is_same_v<Type, Real>) {
    if (!std::isfinite(pivot)) throw SingularMatrixError(...);
}
```
Guards against NaN/Inf propagation — a common source of silent failure.

**Complete solver portfolio:**

| Solver | Specialization | Notes |
|---|---|---|
| GaussJordan | General, produces inverse simultaneously | Full pivoting |
| LU | General (most common) | Partial pivoting, O(n³/3) |
| Cholesky | Symmetric positive definite | 2× faster than LU, checks positive definiteness |
| Band diagonal | Banded systems | O(n·b²) for bandwidth b |
| QR | Ill-conditioned, least squares | `LeastSquaresSolve()` method |
| SVD | Rank-deficient, pseudo-inverse | `Rank()`, `inv_condition()` |
| Iterative | Jacobi, Gauss-Seidel | For large systems with good initial guess |

This covers every solver needed in practice, with appropriate documentation of when to use each.

**Condition number accessibility:**  
`SVDecompositionSolver::inv_condition()` gives the reciprocal condition number directly — a critical diagnostic for ill-conditioned systems that most linear algebra libraries bury in advanced APIs.

---

### 2.5 Field Operations (core/FieldOperations.h)

**Exceptional mathematical documentation:**  
`FieldOperations.h` includes derivation of the gradient formula in spherical coordinates, the physical meaning of divergence and curl, and references to the corresponding textbook sections. This is educational quality documentation inside production library code — rare and valuable.

**Complete coordinate system coverage:**  
`ScalarFieldOperations` supports:
- `GradientCart<N>` — Cartesian
- `GradientSpher` — Spherical (ISO)
- `GradientCyl` — Cylindrical
- `GradientGeneralOrth` — General orthogonal curvilinear (via metric coefficients)

`VectorFieldOperations` provides:
- `DivCart<N>`, `DivSpher`, `DivCyl`
- `CurlCart`, `CurlSpher`, `CurlCyl`
- `JacobianCart<N>`
- `LaplacianCart<N>`, `LaplacianSpher`, `LaplacianCyl`

The general curvilinear version (using metric tensor h₁, h₂, h₃ from `MetricTensor.h`) allows users to define custom coordinate systems and get correct gradient/divergence for free.

---

### 2.6 Coordinate Transformations (core/CoordTransf/)

**ISO 31-11 convention with explicit documentation:**  
`CoordTransfSpherical.h` documents the convention choice (r, θ, φ) with the mathematical definition, warning about the physics vs mathematics convention difference. This avoids the single most common bug in physics code (swapping θ and φ).

**Complete tensor transformation machinery:**  
- Forward transform: `transf(point)` — coordinates
- Inverse: `transfInv(point)` — inverse coordinates  
- Jacobian: `Jacobian(point)` — 3×3 derivative matrix
- Covariant basis vectors: `getBasisVec(point, i)` — ∂r/∂xⁱ
- Unit basis vectors: `getUnitBasisVec(point, i)` — ê_i
- Inverse (contravariant) basis: `getInverseBasisVec(point, i)` — ∇xⁱ
- Covariant vector transform: `transfVecCovariant` — for forces, gradients
- Contravariant vector transform: `transfVecContravariant` — for velocities, displacements

This completeness is graduate-level differential geometry, correctly implemented and accessible via a clean API.

**Numerical pole handling:**  
Division by `sin(θ)` at the poles (θ=0, θ=π) is handled via `clamp` — preventing NaN without silently distorting results far from the poles.

---

### 2.7 AlgorithmTypes.h — Standardized Result Infrastructure

The `AlgorithmStatus` enum with 8 codes (Success, MaxIterationsExceeded, NumericalInstability, SingularMatrix, InvalidInput, Stalled, ToleranceUnachievable, AlgorithmSpecificFailure) is used consistently across the entire library. Combined with `ToString(AlgorithmStatus)`, this creates machine-readable, human-readable failure reporting across all algorithms.

The `ConfigBase` documentation pattern (each algorithm defines its own `struct Config : public ConfigBase`) avoids parameter explosion in function signatures while remaining self-documenting.

---

### 2.8 Richardson Extrapolation and Singularity Handling

`RichardsonExtrapolation.h` provides a standalone extrapolation facility used internally by Romberg integration and optionally by other methods — a correct factoring out of shared infrastructure.

`SingularityHandling.h` provides strategies for near-singular integrals (coordinate splitting, regularization factors, subtraction technique) — addressing one of the most challenging aspects of numerical integration that most libraries leave to the user.

---

### 2.9 Orthogonal Basis Expansions (core/OrthogonalBasis/)

Chebyshev, Hermite, Laguerre, and Legendre bases are classical tools for function approximation, spectral methods, and Gauss quadrature weight computation. Their presence signals MML is targeting serious scientific computing, not just student-level numerics.

---

## 3. Weaknesses

### 3.1 Critical Issues

**W-C1: No automatic differentiation**  
The entire derivative layer is purely numerical. For applications requiring exact gradients (optimization, sensitivity analysis, uncertainty propagation), this is a significant gap. Numerical derivatives:
- Have O(ε^(1/(2k+1))) accuracy for order-2k methods
- Require ≥2 function evaluations per derivative
- Cannot propagate through branching or non-smooth code

Forward-mode automatic differentiation (dual numbers) would provide machine-precision derivatives at the cost of one extra "shadow" float, with zero additional function evaluations.  
**Fix:** Implement `Dual<T>` class for forward-mode AD. This is a medium-sized addition (~300 lines) with large payoff.

---

### 3.2 Significant Issues

**W-S1: Integration `GAUSS10` is the only built-in for Integration2D/3D**  
The higher-dimensional integration routines accept a `GaussRule` enum, but the documentation and examples predominantly show `GAUSS10`. The performance/accuracy tradeoff of different Gauss rules in 2D/3D is not exposed well. For smooth functions in 2D, GAUSS5 may be adequate; for oscillatory functions, more points are needed.  
**Fix:** Add `GAUSS5`, `GAUSS7`, `GAUSS15`, `GAUSS20` labels and document the tradeoff in terms of error order and evaluations.

**W-S2: `IntegrateGaussKronrod` rule selection is opaque**  
The Gauss-Kronrod rules are offered as enum values (`GK_G7K15`, `GK_G10K21`, `GK_G15K31`, etc.) without documentation in the header or AGENTS.md of what to choose. Users who don't already know GK quadrature cannot make an informed choice.  
**Fix:** Add doc comment explaining what each rule trades off (accuracy vs. evaluations) and a recommended default.

**W-S3: `QRSolver` — Householder vs Givens not documented**  
The QR implementation method (Householder reflections vs Givens rotations) affects performance for banded matrices dramatically. If the QR uses Householder (typical for dense matrices), it is suboptimal for banded/sparse cases.  
**Fix:** Document the QR algorithm used, and note that sparse banded systems should use band LU or conjugate gradient instead.

**W-S4: Iterative solver convergence monitoring is limited**  
`LinAlgEqSolvers_iterative.h` provides Jacobi and Gauss-Seidel, which are appropriate for diagonally dominant systems. However, missing are:
- **Conjugate Gradient (CG)** — essential for SPD systems (vastly superior to Jacobi/GS)
- **GMRES** — for non-symmetric systems (generalizes CG)
- **BiCGSTAB** — common alternative to GMRES

**Fix:** Implement at minimum CG and GMRES, which together cover >90% of iterative solver use cases.

**W-S5: `DerivationParametricCurve` — no arc-length parameterization**  
The derivative functions compute all Frenet frame quantities at a given parameter `t`, but the parameter is not necessarily arc-length. For formula involving curvature and torsion in applications (robot path planning, animation), the distinction between the parameterized and arc-length formulations matters significantly.  
**Fix:** Add `reparameterizeByArcLength(a, b, numPoints)` returning a `RealFunction<Real>` that maps arc length `s` to parameter `t`.

---

### 3.3 Minor Issues

**W-M1: `SingularityHandling.h` — weakly documented strategy selection**  
The library provides several singularity handling strategies (subtraction, regularization, coordinate change), but no guide for choosing between them based on the singularity type (integrable or non-integrable, algebraic vs logarithmic).  
**Fix:** Add a doc comment table mapping singularity type to recommended strategy.

**W-M2: `MonteCarloIntegration` — no importance sampling**  
Simple Monte Carlo uses uniform sampling from the integration domain. For functions highly concentrated in small regions, this is inefficient (slow convergence). Importance sampling can reduce variance by orders of magnitude.  
**Fix:** Add `MonteCarloWithImportanceSampling()` variant that accepts a sampling distribution.

**W-M3: `CoordTransf` — no composed transformation chain**  
There is no `ComposedCoordTransf(t1, t2)` that applies `t2 ∘ t1` to go e.g., from cylindrical to spherical. Users must compose manually.  
**Fix:** Implement `ComposedCoordTransf` wrapper that chains `ICoordTransf` objects and combines their Jacobians.

**W-M4: `Curves.h` and `Surfaces.h` naming catalog**  
These files provide named parametric curves (helix, lemniscate, etc.) and surfaces (sphere, torus, etc.). These are excellent for testing and demonstration, but their API is inconsistent — some curves take constructor parameters, others use global factory functions. This inconsistency makes the catalog harder to use programmatically.  
**Fix:** Standardize all curves and surfaces to a constructor-based API with consistent parameter naming.

**W-M5: No `FunctionHelpers` composition operators**  
`FunctionHelpers.h` provides adapters but is missing basic function composition (`compose(f, g)` returning `f∘g`) and binding (`bind(f, x0)` returning `g(t) = f(x0, t)`). These are needed in calculus of variations and optimal control applications.

---

## 4. Architecture Assessment

### 4.1 Layer Coherence

The interfaces → core dependency graph is clean:

```
IFunction.h (pure interfaces — no dependencies)
    ↑
DerivationBase.h (step sizes — depends only on PrecisionValues)
    ↑
Derivation*.h (numerical differentiation — depends on IFunction, DerivationBase)
    ↑
Integration*.h (depends on IFunction, Derivation for adaptive methods)
    ↑
FieldOperations.h (depends on Integration, Derivation, CoordTransf)
    ↑
LinAlgEqSolvers (depends only on Matrix/Vector from base)
```

This is a proper partial order — no cycles, each layer adds capability cleanly.

### 4.2 Standardization Score

The `Config + Result + AlgorithmStatus` pattern is used uniformly in:
- `IntegrationResult` (not exactly Config but has tolerance parameter)
- All LinAlg solvers (implicit config via constructor + explicit solve)
- `GaussKronrod::Result`

Minor inconsistency: the 1D integrators `IntegrateTrap/Simpson/Romberg` accept `tolerance` and `maxIterations` directly as function parameters rather than via a `Config` struct. This slightly inconsistency creates a different call pattern vs. `IntegrateGaussKronrod`.  
**Fix:** Add `IntegrationConfig` struct and offer Config-based overloads.

### 4.3 Missing Module: Spectral / FFT-based integration

For smooth periodic functions, spectral methods (Fourier-based integration) converge super-exponentially. For non-periodic functions, Clenshaw-Curtis quadrature does the same via cosine transforms. Neither is available in the core integration layer (FFT lives in algorithms, and there is no bridge to integration). This is a gap for computational physics applications.

---

## 5. Interface Design Quality

The interface files demonstrate a sophisticated understanding of numerical computing abstractions. The key quality indicators:

1. **Minimal virtual interface** — each interface adds exactly the methods needed for that category, no more.
2. **Return types carry numerical context** — `GetValues()` returns structured (xi, fi) pairs.
3. **Template parameters match mathematical space** — `IScalarFunction<N>` for `ℝᴺ → ℝ`, `IVectorFunctionNM<N,M>` for `ℝᴺ → ℝᴹ`, etc.
4. **Interface segregation** — ODE stepper interface separate from system interface; parametric curve interface separate from function interface.
5. **No god interfaces** — the hierarchy is deep and narrow, not wide.

One concern: `IParametricCurve<N>` inherits `IFunction<VectorN<Real,N>, Real>`. This creates an implicit coupling between the parametric concept and the function interface. For curves with multiple parameterizations (arc-length vs parameter), this coupling prevents easy extension.

---

## 6. Final Grade

| Category | Score (1–10) | Comments |
|---|---|---|
| **Interface richness** | 9.5 | Complete, mathematically coherent, well-segregated |
| **Derivative correctness** | 9.0 | Optimal step sizes, all orders, all function types |
| **Integration completeness** | 9.0 | 1D through 3D, path, surface, improper, Monte Carlo |
| **Linear algebra coverage** | 8.5 | All major solvers; iterative selection could improve |
| **Field operations** | 9.5 | All coordinate systems, excellent documentation |
| **Coordinate transforms** | 9.0 | Full Jacobian, co/contravariant — graduate-level API |
| **Standardization** | 8.0 | Config+Result pattern good; some integrators inconsistent |
| **Missing capabilities** | 7.0 | No AD, no CG/GMRES, no spectral integration |

### Overall Core Layer Grade: **8.7 / 10**

The core layer is the strongest part of MML. The interface design, derivation step size selection, integration breadth, and field operations documentation quality are professional-grade. The main gap is the absence of automatic differentiation and modern iterative solvers (CG, GMRES) — both would significantly expand the library's applicability to large-scale problems.
