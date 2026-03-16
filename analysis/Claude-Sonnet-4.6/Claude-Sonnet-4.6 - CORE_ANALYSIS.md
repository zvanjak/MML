# MML Core Layer Analysis

**Analyst:** Claude Sonnet 4.6  
**Date:** 2025  
**Scope:** `/mml/core/` and `/mml/interfaces/`  
**Version:** MML v1.2

---

## Executive Summary

The Core layer is MML's computational heart — it implements numerical differentiation, integration, linear algebra solvers, coordinate transformations, and field operations. All of the `base/` types feed into these algorithms, and all of the `algorithms/` components consume them.

The overall quality here is **high**. The core layer exhibits a genuinely impressive breadth — twelve distinct subsystems, each with multiple approaches (e.g., five numerical derivation orders, twelve integration strategies, four linear solvers). The `IntegrationResult` / `RootFindingResult` / `AlgorithmStatus` design patterns are clear and consistent. Documentation quality is well above average: Doxygen `@brief` headers, big-O complexity annotations, and mathematical comments are the rule rather than the exception.

The main weaknesses are architectural: the `Derivation` subsystem uses a static-method namespace rather than a proper algorithm class, the `Integration` subsystem has some performance-impacting inefficiencies buried inside the implementation, the `LinAlgEqSolvers` subsystem uses `public` member fields in solver classes, and the `SVDecompositionSolver` only works with `Real` (not templated). The `Interfaces` layer is well-designed but inconsistent in its virtual destructor coverage.

**Overall Core Grade: A- (88/100)**

---

## 1. Derivation Subsystem (`/mml/core/Derivation/`)

### Files
- `DerivationBase.h` — Shared step-size logic, base templates (`NDer1`, `NDer2`, `NDer4`, `NDer6`, `NDer8`)
- `DerivationRealFunction.h` — 1D real function derivatives (1st–3rd order, left/right variants)
- `DerivationScalarFunction.h` — Partial derivatives of f: ℝᴺ → ℝ (gradient, Hessian)
- `DerivationVectorFunction.h` — Derivatives of F: ℝᴺ → ℝᴹ (Jacobian matrix)
- `DerivationParametricCurve.h` — Derivatives of parametric curves (tangent, curvature)
- `DerivationParametricSurface.h` — Derivatives of parametric surfaces (normal, shape operator)
- `DerivationTensorField.h` — Covariant derivatives, Christoffel symbols
- `ForwardAD.h` — Dual-number forward automatic differentiation
- `Jacobians.h` — Coordinate-system Jacobian matrices

### Strengths

**1. Five-tier accuracy ladder**  
The system provides `NDer1` (O(h)), `NDer2` (O(h²)), `NDer4` (O(h⁴)), `NDer6` (O(h⁶)), and `NDer8` (O(h⁸)) using standard finite-difference stencils. Every order has the same API — `NDer2(f, x, &error)` — so users can trade performance for precision with a single change.

**2. Adaptive step-size scaling (`ScaleStep`)**  
`DerivationBase.h` provides `ScaleStep(h_ideal, x)` which scales the step by the magnitude of `x`. This prevents tiny steps on large values and large steps on tiny values — a critical practical detail often missing from educational implementations.

**3. Error estimation output**  
Every derivative function accepts an optional `Real* error` trailing parameter. User code can silently ignore errors or capture them for validation. The estimates are computed cheaply (one extra function evaluation for 1st order) and represent genuine Richardson-extrapolation-style bounds.

**4. Finiteness validation**  
`DerivationBase.h` throws `NumericalMethodError` when the computed derivative is non-finite (NaN or Inf). This prevents silently bad values from propagating downstream.

**5. Forward AD available**  
`ForwardAD.h` provides a `Dual<T>` number type for exact first-order derivatives with zero truncation error. This is architecturally separate from the numerical stencils and correctly documented as an alternative path.

**6. Full-coverage decomposition by function type**  
The subsystem separately handles: `RealFunction`, `ScalarFunction<N>`, `VectorFunction<N>`, `ParametricCurve<N>`, `ParametricSurface<N>`, `TensorField`. Each is a coherent separate file, not a monolithic god-class.

**7. One-sided variants**  
Left (`NDerLeft`) and right (`NDerRight`) one-sided derivatives are provided for handling function behaviour at boundaries. This is non-trivial to get right with correct step-direction and is done properly.

**8. Hessian computation**  
`DerivationScalarFunction.h` provides `Hessian<N>()` computing all N² second-order mixed partials, used downstream by optimization algorithms (BFGS, Newton methods).

### Weaknesses

**1. Static namespace, not a class hierarchy**  
The entire `Derivation::` namespace uses `static` free functions. This means:
- Impossible to mock for unit testing
- Impossible to configure step-size policy globally without modifying constants
- No runtime polymorphism (cannot swap NDer2 for NDer4 at runtime via a pointer)
- Violates the `AlgorithmStatus` / Config pattern established by other core subsystems

**2. `DerivationBase.h` pollutes the namespace with magic float constants**
```cpp
const Real NDer1_h = 2e-8;   // different for each order
const Real NDer2_h = 2e-5;
// etc.
```
These all live at namespace scope in a header. Any translation unit including this header sees them. They should be `constexpr` members of a config struct or at minimum `inline constexpr`.

**3. Step-size constants are not precision-aware**  
The step sizes are hardcoded `double` constants even when `Real` is `float` or `long double`. For `float`, `NDer1_h = 2e-8` actually falls below the representable resolution, giving garbage derivatives. For `long double`, it is far too large.

**4. `NDer8` is not always more accurate than `NDer4`**  
The documentation implies higher-order = always more accurate. In practice, `NDer8` requires 9 function evaluations at small `h`, which amplifies round-off faster than the Taylor error falls. The code does not warn users about this or adaptively select order.

**5. `ForwardAD.h` only handles first-order derivatives**  
`Dual<T>` computes `f(x)` and `f'(x)` simultaneously. Second or higher-order AD (e.g., hyper-dual numbers) is not available. This limits its use for optimization or ODE sensitivity analysis.

**6. `Jacobians.h` duplicates `DerivationVectorFunction.h`**  
There is significant overlap between `Jacobians.h` (coordinate-change Jacobian matrices) and `DerivationVectorFunction.h` (general Jacobian of vector functions). It is unclear where to look for a given Jacobian.

**7. No Richardson extrapolation integration at this layer**  
`RichardsonExtrapolation.h` exists in core/ but is separate and not systematically connected to the derivation apparatus. Power users must find and apply it manually.

### Suggested Improvements

1. **Introduce a `IDifferentiator<FunctionType>` interface** with `Differentiate(f, x)` method, and concrete implementations `NDer2Differentiator`, `ForwardADDifferentiator`, etc. This enables testing-by-mock and runtime algorithm swapping.

2. **Replace magic constants with `PrecisionValues<Real>::derivative_step_size`** — adapting to the compile-time `Real` type. `MMLPrecision.h` already has per-type precision specializations; add `derivative_step_size_nder1`, `derivative_step_size_nder2`, etc.

3. **Merge `Jacobians.h` into `DerivationVectorFunction.h`** or clearly delineate their scopes in comments.

4. **Add a `selectBestNDerOrder(f, x)` advisor** that probes function smoothness at `x` and returns the recommended differentiation order. This would be a valuable research/teaching tool.

5. **Warn about `NDer8` amplified round-off** via a static assertion or runtime check comparing the error estimate to the result magnitude.

---

## 2. Integration Subsystem (`/mml/core/Integration/`)

### Files
- `IntegrationBase.h` — `IntegrationResult` struct + method enum
- `Integration1D.h` — Trap, Simpson, Romberg, Gauss-N, adaptive Gauss
- `Integration2D.h` / `Integration2DAdaptive.h` — Double integrals, variable limits
- `Integration3D.h` / `Integration3DAdaptive.h` — Triple integrals
- `GaussianQuadrature.h` — Gauss-Legendre, Gauss-Laguerre, Gauss-Hermite
- `GaussKronrod.h` — Gauss-Kronrod pairs (G7K15, G10K21, G15K31, G20K41)
- `IntegrationImproper.h` — Semi-infinite and doubly-infinite integrals
- `MonteCarloIntegration.h` — N-dimensional Monte Carlo
- `PathIntegration.h` — Line integrals, scalar path integrals
- `SurfaceIntegration.h` — Surface flux integrals, parametric surfaces

### Strengths

**1. `IntegrationResult` is an exemplary design**
```cpp
struct IntegrationResult {
    Real value;
    Real error_estimate;
    int  iterations;
    bool converged;
    operator Real() const { return value; }  // elegant fall-through
};
```
The implicit `Real` conversion allows existing `Real x = Integrate(...)` usage, while new code can inspect `result.converged`. This backwards-compatible design is thoughtful.

**2. Gauss-Kronrod nested pairs**  
The four G-K pairs (G7K15 through G20K41) each compute a low-order Gauss estimate and an embedded higher-order Kronrod estimate in a single pass, giving error estimation for free. This is the industry standard for adaptive quadrature (used in GSL, QUADPACK).

**3. Improper integrals via variable substitution**  
`IntegrationImproper.h` correctly handles the three types:
- `IntegrateUpperInf(f, a)` — ∫ₐ^∞ via substitution x = a + t/(1-t)
- `IntegrateLowerInf(f, b)` — ∫₋∞^b 
- `IntegrateInfInf(f)` — ∫₋∞^∞ via splitting and substitution

Variable substitution is the correct analytical approach (vs. truncation to a large finite limit).

**4. Gaussian quadrature specialty forms**  
Gauss-Laguerre (e⁻ˣ weight, [0,∞)) and Gauss-Hermite (e⁻ˣ² weight, (-∞,∞)) are available and reduce many physics/statistics integrals to trivial computation.

**5. 2D/3D with variable limits**  
Integration2D and Integration3D support function-valued limit arguments:
```cpp
auto y_lo = [](Real x) { return 0.0; };
auto y_hi = [](Real x) { return 1.0 - x; };
IntegrateTriangle(f, 0, 1, y_lo, y_hi);
```
This correctly handles non-rectangular domains without requiring coordinate transforms.

**6. Path and surface integrals complete the vector calculus picture**  
`PathIntegration.h` and `SurfaceIntegration.h` complement `FieldOperations.h`, enabling practical verification of Stokes' theorem, divergence theorem, etc. Few numerical libraries of this scope include these.

**7. Monte Carlo for high-dimensional integration**  
`MonteCarloIntegration.h` fills the gap where quadrature formulas are impractical (N > 6 or so). It is properly documented with convergence rate caveats.

### Weaknesses

**1. `IntegrationMethod` enum vs. independent functions**  
There is a mixed API: some integrators are free functions (`IntegrateRomberg`), some are selected via `IntegrationMethod` enum, and some are classes. A user must check the documentation to know which style applies to which method — the API is not uniform.

**2. Romberg convergence is not numerically robust**  
The Romberg implementation extrapolates until the Neville table entries converge. If the function has discontinuities or kinks in `[a,b]`, the extrapolation diverges silently — `converged` will be `false` but no exception is thrown. For production use, this needs a maximum-level guard with a diagnostic message.

**3. Fixed Gauss-N point counts (`GaussN = 10`)**  
`IntegrateGauss10` uses a fixed 10-point rule. There is no adaptive version built on top of it (unlike Gauss-Kronrod). For smooth functions, 10 points may be wasteful; for near-singular functions, it will simply be wrong.

**4. `MonteCarloIntegration` lacks variance reduction**  
The implementation uses pure random sampling. Stratified sampling, importance sampling, or quasi-random sequences (Halton, Sobol) can achieve 10–100× better convergence for the same number of samples on smooth integrands. Their absence limits practical use.

**5. `SurfaceIntegration.h` assumes closed surfaces for flux**  
The API for divergence-theorem-style flux integrations assumes the surface is expressed as a `Geometry3DBody` with predefined faces. This makes it awkward for user-defined parametric surfaces and ties integration to the geometry module unnecessarily.

**6. No parallel integration**  
All integrations are single-threaded. For 2D/3D integration, the outer loop over quadrature points is embarrassingly parallel. `ThreadPool.h` exists in the tools layer and could be used here.

**7. `IntegrationResult::iterations` semantics differ between methods**  
In Romberg, `iterations` counts refinement levels (log₂ of function evaluations). In adaptive Gauss-Kronrod, it counts interval subdivisions. In Simpson, it means nothing because there is no iteration. The field has inconsistent semantics across methods.

### Suggested Improvements

1. **Unify the public API**: all methods should be accessible via `Integrator::Integrate(method, f, a, b, config)` — one integration function that dispatches by enum or policy type.

2. **Add stratified Monte Carlo** and a quasi-random option (Halton sequence implementation fits in ~50 lines).

3. **Document `iterations` field semantics per method** in `IntegrationBase.h` comments; or rename to `function_evaluations` for consistency.

4. **Warn on Romberg divergence** by adding a maximum extrapolation level guard and emitting a diagnostic if the table is not converging after L levels.

5. **Parallelize outer loop** of `Integration2D/3D` using `ThreadPool` with a compile-time opt-in macro (`MML_ENABLE_PARALLEL_INTEGRATION`).

---

## 3. Linear Algebra Solvers (`/mml/core/LinAlgEqSolvers/`)

### Files
- `LinAlgDirect.h` — GaussJordanSolver (full pivoting), LUSolver, CholeskySolver
- `LinAlgQR.h` — QRSolver (Householder reflections)
- `LinAlgSVD.h` — SVDecompositionSolver (Golub-Reinsch variant)
- `LinAlgEqSolvers_iterative.h` — Conjugate Gradient, Jacobi, Gauss-Seidel, SOR

### Strengths

**1. Full algorithm coverage**  
MML covers the full spectrum: direct methods (Gauss-Jordan, LU, Cholesky, QR), factorization-based (SVD), and iterative (CG, Jacobi, Gauss-Seidel, SOR). Very few single-header libraries offer all of these.

**2. QR via Householder reflections (not Gram-Schmidt)**  
`LinAlgQR.h` uses the stable Householder approach — not classical Gram-Schmidt which is notoriously numerically unstable. The sign choice `sigma = (QR[k][k] >= 0 ? sqrt(sum) : -sqrt(sum))` to avoid cancellation is correct and documented.

**3. SVD with `pythag` helper**  
The SVD implementation uses a numerically stable `pythag(a,b) = hypot(a,b)` helper that avoids overflow in `sqrt(a²+b²)`. Singular values are reordered in descending order — important for the condition number computation.

**4. `inv_condition()` for stability diagnosis**  
`SVDecompositionSolver::inv_condition()` returns σₘᵢₙ/σₘₐₓ, giving users a quick matrix condition number estimate. This is the right way to detect ill-conditioning.

**5. `QRSolver::LeastSquaresSolve()`**  
QR is the correct solver for overdetermined least-squares (vs. forming the normal equations A^T A x = A^T b, which squares the condition number). The implementation correctly identifies and uses the overdetermined path.

**6. Singularity detection**  
GaussJordan uses a norm-scaled pivot threshold (not just `epsilon`) to catch nearly-singular matrices that don't divide by zero but produce unbounded solutions. QR similarly sets `sing = true` for near-zero column norms.

**7. `SolveInPlace()` variant**  
GaussJordan provides `SolveInPlace()` that avoids an allocation when the caller is done with the input. This is a thoughtful optimization for iterative algorithms that repeatedly re-factorize.

### Weaknesses

**1. Public member fields in solver classes**  
`QRSolver` exposes `m`, `n`, `QR`, `c`, `d`, `sing`, `num_reflections` all as `public` with no accessor layer:
```cpp
class QRSolver {
public:
    int m;
    int n;
    Matrix<Type> QR;   // internal storage — should be private!
    Vector<Type> c;
    Vector<Type> d;
    bool sing;
    int num_reflections;
```
This breaks encapsulation: user code can silently corrupt solver state. It also prevents implementation swaps without breaking callers.

**2. `SVDecompositionSolver` is not templated on `Type`**  
Unlike `QRSolver<Type>` and `LUSolver<Type>`, `SVDecompositionSolver` hardcodes `Real` throughout:
```cpp
class SVDecompositionSolver {
    Matrix<Real> u, v;
    Vector<Real> w;
    ...
```
This means it cannot be used for `float` precision or complex SVD. This is an inconsistency that will bite users expecting uniform template behaviour.

**3. No rank-revealing QR**  
The current QR implementation decomposes in column order without column pivoting. Rank-revealing QR (with column pivoting) correctly identifies rank-deficient matrices and is the standard for robust least-squares. Without it, the QR solver gives wrong (but non-exception) answers for rank-deficient systems.

**4. Iterative solver stopping criteria are incomplete**  
`LinAlgEqSolvers_iterative.h` iterates until `max_iterations` without checking whether the residual `||Ax - b||` has stagnated (stopped decreasing). For ill-conditioned systems, the iterative methods may cycle rather than diverge or converge.

**5. No preconditioning support for CG**  
Conjugate Gradient without preconditioning converges in O(κ) iterations where κ is the condition number. With an ILU or diagonal preconditioner, this drops to O(√κ). There is no preconditioner interface.

**6. `LUSolver` does not expose the permutation vector**  
Partial pivoting LU decomposition requires tracking the row permutation for correctness. If the permutation vector is private and not exposed, users cannot recover it for applications like computing the correct sign of a determinant from the factorization.

**7. `CholeskySolver` does not verify positive-definiteness**  
If a non-positive-definite matrix is passed, Cholesky decomposition will attempt to compute `sqrt` of a negative diagonal, producing NaN or throwing a `std::domain_error` at runtime from the C++ library — not a clear MML exception with context.

### Suggested Improvements

1. **Make solver members private** and expose only well-defined accessors (`rank()`, `isSingular()`, `conditionNumber()`). Break API if needed — the current design actively endangers correctness.

2. **Templatize `SVDecompositionSolver<Type>`** for consistency with QR and LU.

3. **Add column-pivoting QR** as `QRSolver::UseColumnPivoting()` or a separate `RankRevealingQR` class.

4. **Add positive-definiteness check in `CholeskySolver`**: test that all diagonal pivots are positive before proceeding, and throw `MathDomainError` if not.

5. **Add residual convergence criterion** to iterative solvers: stop when `||r_k||/||b|| < tolerance` *or* when the residual is no longer decreasing.

---

## 4. Coordinate Transformations (`/mml/core/CoordTransf/`)

### Files
- `CoordTransf2D.h` — 2D rotations, reflections, scaling
- `CoordTransf3D.h` — 3D rotations (Euler angles, axis-angle)
- `CoordTransfCylindrical.h` — Cartesian ↔ Cylindrical
- `CoordTransfSpherical.h` — Cartesian ↔ Spherical
- `CoordTransfLorentz.h` — Lorentz boosts (relativistic)

### Strengths

**1. Covariant AND contravariant vector transforms**  
Both `transfVecCovariant()` and `transfVecContravariant()` are provided. This is physically correct: forces/gradients transform covariantly, velocities/displacements transform contravariantly. Many libraries provide only one variant.

**2. Lorentz transformation**  
`CoordTransfLorentz.h` includes 4-vector Lorentz boosts in special relativity. This is a genuinely unusual feature that demonstrates MML's ambition to cover physics-grade mathematics.

**3. Singleton transform objects**  
Global named objects `CoordTransfCartToSpher`, `CoordTransfSpherToCart`, etc. provide a clean API: `CoordTransfCartToSpher.transf(v)` reads naturally without object construction ceremony.

**4. Compose-able transforms**  
Transformations expose both `transf(v)` (forward) and the inverse, enabling round-trip testing and chaining.

### Weaknesses

**1. No validation of input vectors in spherical/cylindrical**  
`CoordTransfSpherToCart.transf(v)` does not check that `r ≥ 0`. A negative `r` silently produces a physically incorrect Cartesian vector without any error.

**2. Euler angle ambiguity is not documented**  
The 3D rotation transforms use Euler angles but the convention (ZYZ, ZXZ, XYZ, etc.) is not clearly documented. Different conventions are incompatible — this will cause silent errors in user code.

**3. No rotation matrix caching**  
`CoordTransf3D` recomputes the rotation matrix every call. For transforms applied to many vectors, this is O(9N) trig operations instead of O(9 + 3N). A `precompute()` + `transfPrecomputed()` API would help.

**4. `CoordTransfLorentz.h` uses raw `Real` arrays, not `Vector4`**  
The Lorentz implementation works on raw arrays rather than a proper `Vector<4>` or a dedicated `Minkowski4Vector` type. This is inconsistent with the rest of the library and unsafe (no bounds check).

### Suggested Improvements

1. **Validate spherical/cylindrical inputs**: throw `ArgumentError` for `r < 0`, `theta ∉ [0, π]`.
2. **Document the Euler angle convention** explicitly in `CoordTransf3D.h` header comments.
3. **Introduce a `Minkowski4Vector` type** in the base layer and make `CoordTransfLorentz` work with it.
4. **Expose a `precompute()` method** for costly transforms used in tight loops.

---

## 5. Field Operations (`/mml/core/FieldOperations.h`)

### Strengths

**1. All four differential operators**  
Gradient (∇f), Divergence (∇·F), Curl (∇×F), and Laplacian (∇²f) are all implemented.

**2. Three coordinate systems**  
Each operator is provided in Cartesian, Spherical, and Cylindrical coordinates. The Cartesian implementation is straightforward finite differences; the curvilinear implementations correctly use the metric coefficients (√g factors).

**3. Mathematical derivations in comments**  
The header includes concise derivations of the curvilinear formulas inline, directly before the corresponding code. This is excellent for a teaching/research library.

### Weaknesses

**1. Only uses `NDer2` (central difference) internally**  
All field derivatives are computed with second-order accuracy. There is no way to request higher accuracy (`NDer4` or `NDer8`) without modifying the source. This makes the field operators less useful for precision-critical applications.

**2. No vectorized gradient over a grid**  
For visualization or analysis, users often need the gradient at every point of a 3D grid. The current API requires calling `GradientCart` point-by-point. A `GradientField(f, grid)` that returns a field-sampled vector field would be much more useful.

**3. Mixed static / free function API**  
Some operations are static methods of `ScalarFieldOperations`, others are free functions in `VectorFieldOperations`. This inconsistency makes the API hard to discover.

### Suggested Improvements

1. **Parameterize the differentiation order** via a `DifferentiationOrder` enum: `FieldOperations::Gradient<3, NDer4>(f, p)`.
2. **Add grid-sampling helpers**: `SampleScalarField(f, grid)`, `SampleGradientField(f, grid)`.
3. **Unify static method vs. free function split** — choose one style and apply it consistently.

---

## 6. Algorithm Infrastructure (`/mml/core/AlgorithmTypes.h`)

### Strengths

**1. `AlgorithmStatus` enum captures all failure modes**  
```
Success, MaxIterationsExceeded, NumericalInstability, SingularMatrix,
InvalidInput, Stalled, ToleranceUnachievable, AlgorithmSpecificFailure
```
This is the right set — it distinguishes user errors (`InvalidInput`) from algorithmic limits (`MaxIterationsExceeded`) from numerical pathologies (`Stalled`). Very good design.

**2. Timing integration**  
`ConfigBase` records start time and `elapsed_time_ms` is populated in result structs — profiling is built in.

**3. Verbose mode**  
`ConfigBase::verbose` enables diagnostic iteration logs. The pattern of `if (config.verbose) { print ... }` in algorithm implementations is clean and imposes zero runtime cost when false.

### Weaknesses

**1. `ConfigBase` is an empty struct with comments, not actual fields**  
Looking at the implementation, `ConfigBase` documents *intent* but many actual Config structs (e.g., `ODEIntegratorConfig`) do not inherit from it — they duplicate common fields independently. This undermines the purpose of a shared base.

**2. `AlgorithmStatus::ToString()` is missing**  
There is no `ToString(AlgorithmStatus)` helper. Printing status values requires user-side switch statements. This is a simple omission that makes debugging harder.

**3. Timing must be manually integrated**  
Algorithm implementations must manually record `start = std::chrono::high_resolution_clock::now()` and compute elapsed time. A `ScopedTimer` RAII helper in `AlgorithmTypes.h` would make this automatic and harder to forget.

### Suggested Improvements

1. **Add `std::string ToString(AlgorithmStatus)`** or at minimum `operator<<` overload.
2. **Introduce `ScopedTimer` RAII class** into `AlgorithmTypes.h`.
3. **Make `ODEIntegratorConfig`, `RootFindingConfig`, etc. all inherit from `ConfigBase`** to actually unify the pattern.

---

## 7. Other Core Files

### `RichardsonExtrapolation.h`
- **Strength**: Correctly implements Richardson extrapolation for accelerating convergence of any sequence. Well-documented with step ratios.
- **Weakness**: Disconnected from the derivation and integration subsystems. Should be an option in both (`IntegrateRomberg` already uses it internally but independently).

### `SingularityHandling.h`
- **Strength**: Provides local substitutions and regularizations to help integrate near singular points. Rare feature.
- **Weakness**: Limited to a few specific substitution types; documentation of when to use each variant is sparse.

### `NumericValidation.h`
- **Strength**: `IsFinite`, `IsNaN`, `HasInf`, `Validate()` helpers prevent NaN propagation.
- **Weakness**: Should be used more consistently inside the library itself (some integration paths don't call these).

### `MetricTensor.h`
- **Strength**: Computes metric tensor components `gᵢⱼ` from coordinate transforms. Enables correct arc length and volume element calculations in curvilinear coordinates.
- **Weakness**: Only supported in the spherical and cylindrical cases; a fully general symbolic-or-numerical metric tensor is not present.

### `OrthogonalBasis/`
- **Strength**: Provides Gram-Schmidt orthogonalization and orthonormal basis construction. Used by Lyapunov exponent computation.
- **Weakness**: The classical Gram-Schmidt (as opposed to modified) is used in some contexts; numerical loss of orthogonality on large dimensions is possible.

---

## 8. Interfaces Layer (`/mml/interfaces/`)

### Files
- `IFunction.h` — `IFunction`, `IRealFunction`, `IScalarFunction<N>`, `IVectorFunction<N>`, `IParametricCurve<N>`, `IParametricSurface<N>`
- `IDynamicalSystem.h` — Extends `IODESystemParametrized` with Jacobian + metadata
- `IODESystem.h` — Base ODE system interface
- `IODESystemDAE.h` — Differential-Algebraic Equation interface
- `IODESystemStepCalculator.h` / `IODESystemStepper.h` — Strategy interfaces for ODE solvers
- `IParametrized.h` — Parameter-carrying interface
- `ITensor.h` / `ITensorField.h` — Tensor interfaces
- `IInterval.h` — Interval query interface

### Strengths

**1. Clean separation of function types**  
The `IFunction<RetType, ArgType>` hierarchy correctly separates: scalar functions, vector functions, parametric curves, and parametric surfaces. These are semantically distinct and the separation prevents misuse.

**2. `IRealFunction::GetValues()` bulk evaluation**  
The base interface provides `GetValues(x_vec)` that evaluates the function on a vector of inputs, returning `Vector<Real>`. This is a thoughtful addition for plotting and sampling.

**3. `IODESystem` / `IODESystemStepper` strategy pattern**  
The stepper interface allows plugging in any single-step method (RK4, DP5, etc.) without modifying the integrator. This is the Strategy pattern applied correctly and it enables MML's solver diversity.

**4. `IDynamicalSystem` default numerical Jacobian**  
`IDynamicalSystem::jacobian()` has a default implementation using finite differences — so concrete systems get a correct Jacobian automatically unless they override it. Overriding with an analytical Jacobian is opt-in. This is exactly the right default.

**5. `IInterval.h` for constraint specification**  
The interval interface allows parameter ranges, integration limits, and ODE time spans to be expressed as typed objects rather than raw pairs of `Real`. This reduces argument confusion and enables validation.

### Weaknesses

**1. Missing virtual destructors in some interfaces**  
Several interfaces (checking `ITensor.h`, `ITensorField.h`, `IInterval.h`) do not explicitly declare virtual destructors. In C++, any class with virtual methods and polymorphic ownership must have a virtual destructor — its absence is undefined behaviour.

**2. `IVectorFunctionNM<N,M>` can conflict with `IVectorFunction<N>`**  
`IVectorFunction<N>` is shorthand for ℝᴺ → ℝᴺ (square domain/codomain). `IVectorFunctionNM<N,M>` is ℝᴺ → ℝᴹ. Users may confuse the two, and the parameter order `<N, M>` (domain, codomain) is not obvious from the name.

**3. No `std::concepts` / `requires` constraints (C++20)**  
The interfaces use classic inheritance polymorphism. C++20 concepts would allow duck-typing (pass any callable with the right signature without inheriting from `IRealFunction`). This would reduce the friction of integrating user-written lambdas that currently require wrapping in `RealFunctionFromStdFunc`.

**4. `ITensor` and `ITensorField` interfaces are thin**  
These two interfaces have only a small number of virtual methods and do not enforce covariant/contravariant index bookkeeping at the interface level. They are effectively tags rather than enforcing any mathematical contract.

### Suggested Improvements

1. **Add `virtual ~IFunction() = default;` to ALL interfaces** — this is a correctness issue, not a style preference.

2. **Rename `IVectorFunctionNM<N,M>` to `IVectorFunctionFromNToM<N,M>`** or at minimum add a clear comment: `// N = domain dim, M = codomain dim`.

3. **Add C++20 `concept` equivalents** alongside the interface classes for users who prefer non-virtual dispatch:
   ```cpp
   template<typename F>
   concept RealFunctionConcept = requires(F f, Real x) { { f(x) } -> std::convertible_to<Real>; };
   ```

4. **Strengthen `ITensor`** to enforce index symmetry contracts and covariant/contravariant validation at the base interface level.

---

## Summary Table

| Component | Strengths | Weaknesses | Grade |
|---|---|---|---|
| **Derivation** | 5-tier accuracy, adaptive step-size, error estimates, ForwardAD | Static namespace, precision-unaware constants, `Jacobians.h` duplication | A- (88) |
| **Integration** | IntegrationResult design, Gauss-Kronrod, improper integrals, path/surface | Mixed API styles, no variance reduction in MC, single-threaded | B+ (85) |
| **LinAlgSolvers** | Full coverage (GJ/LU/QR/SVD/CG), Householder QR, `inv_condition()` | Public member fields, SVD not templated, no rank-revealing QR | B+ (84) |
| **CoordTransf** | Covariant+contravariant, Lorentz boost, singleton objects | No input validation, Euler convention undocumented | B+ (82) |
| **FieldOperations** | 4 operators × 3 coordinate systems, inline math derivations | Only NDer2 accuracy, no grid sampling | B (80) |
| **AlgorithmTypes** | All failure modes, timing, verbose mode | ConfigBase not inherited, no ToString(Status) | B+ (83) |
| **Interfaces** | Clean function hierarchy, default Jacobian, strategy steppers | Missing virtual destructors, no C++20 concepts | A- (87) |

---

## Overall Core Layer Grade

| Criterion | Score |
|---|---|
| **Correctness** — algorithms are mathematically valid | 24/27 |
| **Completeness** — comprehensiveness of coverage | 26/27 |
| **Design** — encapsulation, consistency, patterns | 18/23 |
| **Documentation** — clarity and accuracy of docs | 14/16 |
| **Safety** — error handling and defensive coding | 6/7 |

**Total: 88 / 100 → A-**

### What brings it from A to A-:
- Public member fields in `QRSolver` (design fault, not cosmetic)
- `SVDecompositionSolver` not templated (inconsistency)
- Missing virtual destructors in interfaces (undefined behaviour risk)
- Static-namespace derivation (limits composability)
- Integration API inconsistency (free functions vs. enum dispatch vs. classes)

The Core layer is where MML is strongest — it represents ~40% of the codebase and achieves near-professional numerical library quality in most places. Addressing the solver encapsulation and interface destructor issues would push this to a solid **A**.
