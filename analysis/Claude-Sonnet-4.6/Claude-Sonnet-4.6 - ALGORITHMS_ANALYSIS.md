# MML Algorithms & Systems Analysis

**Analyzed by:** Claude Sonnet 4.6  
**Date:** 2026-03-13  
**Scope:** `mml/algorithms/` (all subdirectories) and `mml/systems/` (all files)

---

## 1. Scope Overview

The algorithms and systems layers sit above the core mathematical infrastructure and provide complete, ready-to-use numerical algorithms and dynamical system analysis tools:

| File / Folder | Purpose |
|---|---|
| `algorithms/ODESolvers/` | Adaptive, fixed-step, stiff ODE integrators; BVP shooting |
| `algorithms/DAESolvers/` | Backward Euler, BDF2/4, RadauIIA, RODAS for DAE systems |
| `algorithms/RootFinding/` | Bisection, Brent, Newton, Secant, Ridders; polynomial roots |
| `algorithms/Optimization.h` | 1D minimization: GoldenSection, Brent, BracketMinimum |
| `algorithms/OptimizationMultidim.h` | N-D: NelderMead, Powell, CG, BFGS, Levenberg-Marquardt |
| `algorithms/EigenSystemSolvers.h` | Jacobi (symmetric), QR iteration (general), Power method |
| `algorithms/Statistics.h` | Descriptive stats: mean, variance, moments, median, t-tests |
| `algorithms/Fourier.h` | DFT (O(n²)), FFT (Cooley-Tukey), DCT, DST |
| `algorithms/Distributions.h` | Probability distributions (PDF, CDF, inverse) |
| `algorithms/CurveFitting.h` | Least squares, polynomial and nonlinear fitting |
| `algorithms/ChebyshevApproximation.h` | Minimax polynomial approximation |
| `algorithms/CompGeometry/` | Convex hull, KD-tree, Voronoi, Delaunay triangulation |
| `algorithms/FieldAnalyzers.h` | Field visualization analysis |
| `algorithms/FieldLineTracer.h` | Streamlines and field line computation |
| `algorithms/FunctionsAnalyzer.h` | Zero/extrema/inflection point detection, asymptotics |
| `algorithms/MatrixAlg.h` | Matrix polynomial, functions (matrix exp, log, sqrt) |
| `algorithms/GraphAlgorithms.h` | Shortest path, spanning tree, flow algorithms |
| `algorithms/GraphSpectral.h` | Graph Laplacian, spectral clustering |
| `systems/DynamicalSystemBase.h` | Template base class for all N-dimensional dynamical systems |
| `systems/DynamicalSystemAnalyzers.h` | Fixed point finder, Lyapunov exponents, bifurcation analysis |
| `systems/ContinuousSystems.h` | Lorenz, Rössler, Van der Pol, Duffing, Hodgkin-Huxley, etc. |
| `systems/DiscreteMaps.h` | Logistic map, Hénon map, standard map |
| `systems/LinearSystem.h` | Linear dynamical systems |
| `systems/DynamicalSystemTypes.h` | Result types: FixedPoint, LyapunovResult, BifurcationResult |

---

## 2. Strengths

### 2.1 Adaptive ODE Integrator

**FSAL (First Same As Last) optimization:**  
`ODEAdaptiveIntegrator<DormandPrince5_Stepper>` implements the FSAL property of the Dormand-Prince method: the last function evaluation of one step is reused as the first of the next. For a 5-stage method, this effectively gives a 4-stage cost per accepted step — a significant performance gain.

**Dense output via Hermite interpolation:**  
When `output_interval < step_size`, the integrator outputs intermediate values by Hermite polynomial interpolation between accepted steps. This is the mathematically correct approach to "continuous" ODE output, vital for event detection and visualization at arbitrary resolution.

**Comprehensive `SolutionStatistics`:**  
```cpp
struct SolutionStatistics {
    int accepted_steps;
    int rejected_steps;
    int function_evaluations;
    Real min_step_size;
    Real max_step_size;
    Real elapsed_time_ms;
    AlgorithmStatus status;
};
```
Users can quantify step size adaptation behavior, identify stiffness (high rejection rate), and profile function evaluation cost — all from one result struct.

**`ODEIntegratorConfig` factory methods:**  
```cpp
static ODEIntegratorConfig Stiff()    // tight tolerances, reduced max_step
static ODEIntegratorConfig Fast()     // relaxed tolerances, larger steps
```
These "preset" factories follow the principle of "provide sensible defaults for common use cases." Users solving standard non-stiff systems don't need to tune tolerances; stiff systems get a safe starting point.

**Step size adapter correctness:**  
The adaptive step control uses `h_new = h * min(max_scale_factor, max(min_scale_factor, safety_factor * (tol / error)^(1/p)))` with the theoretically correct exponent `1/(p+1)` for order-p methods. The safety factor (typically 0.9) and scale limits prevent pathological step size oscillations.

**Full stepper library:**

| Stepper | Type | Best For |
|---|---|---|
| `ODESystemStepperEuler` | Fixed, order-1 | Educational |
| `ODESystemStepperRK4` | Fixed, order-4 | Simple non-stiff problems |
| `DormandPrince5_Stepper` | Adaptive, order-5 | General non-stiff |
| `DormandPrince8_Stepper` | Adaptive, order-8 | High accuracy |
| `CashKarp_Stepper` | Adaptive, order-5 | Alternative to DP5 |
| `BulirschStoer_Stepper` | Adaptive, extrapolation | Very smooth solutions |

This is a complete portfolio covering the spectrum from simple teaching cases to high-accuracy production computation.

---

### 2.2 DAE Solvers — A Rare Achievement

**High-index DAE support via multiple methods:**  
The presence of `DAEBDF2`, `DAEBDF4`, `DAERadauIIA`, and `DAERODAS` in a single-header library is remarkable. Most open-source numerical libraries either have no DAE support or provide only index-1 backward Euler:

| Solver | Method | Order | Suitable For |
|---|---|---|---|
| `DAEBackwardEuler` | Implicit Euler | 1 | Index-1 DAEs, stiff, first iteration |
| `DAEBDF2` | Backward differentiation formula | 2 | Mildly stiff DAEs |
| `DAEBDF4` | BDF order 4 | 4 | Moderately stiff, good accuracy |
| `DAERadauIIA` | Radau IIA implicit Runge-Kutta | 5 | Stiff, high accuracy |
| `DAERODAS` | Rosenbrock-Wanner | 4 | Stiff, "black box" robustness |

BDF4 is the industry standard in circuit simulation (SPICE) and chemical reactor modeling. RadauIIA is the solver of choice for stiff problems requiring high accuracy. Having both alongside a robust Rosenbrock method gives MML coverage that parallels MATLAB's `ode15s`.

---

### 2.3 Root Finding

**`FindRootBrent` — provably convergent, superlinear:**  
Brent's method combines the reliability of bisection (guaranteed convergence) with the speed of inverse quadratic interpolation / secant method when possible. The implementation correctly falls back to bisection when the faster step would violate the bracket, making it impossible to diverge.

**Newton-Raphson with numerical gradient and bracket fallback:**  
The Newton implementation uses `NDer4` for the numerical derivative, automatically correct without user-provided gradients. When Newton steps escape the bracket `[x1, x2]`, it falls back to bisection — preventing the classic Newton divergence scenario on non-convex functions. This is a practical improvement over pure Newton.

**Dual API — ergonomic design:**  
```cpp
// Simple API
Real root = RootFinding::FindRootBrent(f, a, b, 1e-10);

// Full API
RootFindingResult result = RootFinding::FindRootBrent(f, a, b, config);
```
The simple API is clean for interactive use and educational code. The full API gives convergence history, iteration count, and AlgorithmStatus for production code. This dual pattern is consistently applied across all root-finding methods.

**`FindRootBrackets` — multi-root discovery:**  
The function sweeps an interval for sign changes, returning multiple `[lo, hi]` bracket pairs. This is the correct first step for multi-root problems and is often missing from single-root-focused libraries.

**Polynomial root finding:**  
`SolveQuadratic`, `SolveCubic`, `SolveQuartic` return complex roots and handle all degenerate cases (discriminant = 0, complex conjugate pairs). The numeric closures handle `b²−4ac → 0` edge cases without cancellation errors.

---

### 2.4 Optimization

**`BrentMinimize` — gold standard for 1D:**  
The golden-section + parabolic interpolation combination correctly handles flat regions (where parabolic interpolation would divide by zero) by falling back to golden section. Implementation includes the "acceptance criterion" ensuring the step doesn't backtrack.

**`NelderMead` — strong implementation:**  
The simplex-based minimizer correctly implements reflection, expansion, contraction, and shrinkage operations with their standard coefficients (α=1, γ=2, ρ=0.5, σ=0.5). The `MultidimMinimizationResult` includes timing and status.

**`LevenbergMarquardt` — for nonlinear least squares:**  
The Levenberg-Marquardt algorithm fuses gradient descent and Gauss-Newton, making it the standard for nonlinear least squares (curve fitting, model parameter estimation). Its inclusion distinguishes MML from pure optimization libraries that focus only on unconstrained minimization.

**`BFGS` — quasi-Newton for smooth functions:**  
BFGS approximates the Hessian from gradient differences, achieving superlinear convergence for smooth objectives without second-order derivatives. The presence of BFGS alongside NelderMead covers both smooth and non-smooth landscapes.

---

### 2.5 Eigenvalue Solvers

**`SymmMatEigenSolverJacobi` — guaranteed convergence:**  
The Jacobi algorithm for symmetric matrices is one of the most reliable eigensolvers — it convergences unconditionally (quadratic rate) and produces a full orthonormal eigenvector matrix. For small-to-medium symmetric matrices (physics, mechanics, chemistry), this is the correct choice.

**Comprehensive `EigenSolverConfig`:**  
The config struct provides `HighPrecision()` and `Fast()` factory presets alongside individual knobs (`tolerance`, `max_iterations`, `sort_eigenvalues`, `compute_eigenvectors`, `verbose`). The `sort_eigenvalues` flag is practical — many applications assume ascending order.

**`EigenSolverResult` quality:**  
Includes `residual` (‖AV − VΛ‖), `elapsed_time_ms`, `AlgorithmStatus` — giving users a complete picture of solution quality.

---

### 2.6 Statistics

**Two-pass stable variance:**  
`AvgVar(data)` uses the two-pass formula: first compute mean, then sum squared deviations from mean. This avoids the catastrophic cancellation of the naive single-pass `Σx² − (Σx)²/n` formula that fails when variance ≪ mean². This is a correctness feature that separates professionally written numerical code from naive implementations.

**`Moments` struct:**  
Returns mean, average absolute deviation, standard deviation, variance, skewness, and kurtosis in a single pass. Skewness and kurtosis require special normalization (corrected for bias in the kurtosis denominator). These are correct.

**`Percentile`:**  
Uses linear interpolation between sorted values rather than the nearest-rank method — the standard for continuous distributions.

---

### 2.7 Fourier Analysis

**Dual FFT/DFT offering:**  
`DFT` (O(n²), any size) alongside `FFT` (Cooley-Tukey, power-of-2 only) is the correct combination: use DFT for verification or odd sizes, FFT for production performance.

**`FourierValidation` namespace:**  
A dedicated namespace providing input validators (dimension check, power-of-2 check) shared between DFT, FFT, DCT, DST. This factoring-of-validation prevents the common pattern of copy-pasted checks across variants.

**DCT-II and DCT-III (forward and inverse):**  
These are the variants used in JPEG compression, audio coding (MP3, AAC), and spectral methods for smooth functions. Their presence alongside the standard FFT signals that MML targets signal processing applications beyond basic spectrum analysis.

---

### 2.8 Computational Geometry (algorithms/CompGeometry/)

**Comprehensive 2D and 3D coverage:**  
- Convex hull (2D), `ConvexHull3D`
- Delaunay triangulation, Voronoi diagram
- KD-tree (point proximity, range queries)
- Polygon operations, intersection tests

`KDTree` enables O(log n) nearest-neighbor queries — essential for large-scale geometry, mesh generation, and spatial hashing. Its presence in MML is unusually complete for a mathematical library.

---

### 2.9 Dynamical Systems Layer (mml/systems/)

**`DynamicalSystemBase<N,P>` — elegant template design:**  
The template encodes state dimension `N` and parameter count `P` at compile time:
```cpp
class LorenzSystem : public DynamicalSystemBase<3, 3> { ... };
```
State names, parameter names, and parameter ranges are stored at runtime — enabling generic analysis tools to operate on any specific system without case-specific knowledge.

**`FixedPointFinder` — Newton + LU:**  
Uses Newton-Raphson (with numerically computed Jacobian via `DerivationVectorFunction`) and LU decomposition to find fixed points. `FindMultiple()` accepts a list of initial guesses and filters duplicates (within tolerance). The uniqueness filtering prevents reporting the same fixed point twice from closely clustered initial conditions — a practical detail.

**`ClassifyFixedPoint()` — eigenvalue-based classification:**  
Fixed point stability is determined correctly from eigenvalues of the Jacobian at the fixed point:
- All eigenvalues Re < 0 → stable node/spiral
- At least one Re > 0 → unstable
- Purely imaginary → center
- Mixed real/imaginary → saddle
- Complex conjugate pairs → spiral

This is the rigorous mathematical criterion (Hartman-Grobman linearization theorem). The implementation also labels the type (Node, Spiral, Saddle, Center, StarNode) — providing physically meaningful vocabulary.

**`LyapunovAnalyzer` — Gram-Schmidt orthonormalization:**  
Lyapunov exponents are computed via periodic Gram-Schmidt orthonormalization of the perturbation vectors (the standard Benettin algorithm). The result includes the full spectrum plus:
- `maxExponent` — positive implies chaos
- `kaplanYorkeDimension` — fractal dimension from the Kaplan-Yorke formula

This is graduate-level dynamical systems analysis correctly implemented.

**`BifurcationAnalyzer::Sweep()` — full bifurcation diagrams:**  
The analyzer sweeps a parameter over a range, runs each trajectory to steady state, records Poincaré section maxima, and returns the full diagram data. This is the standard computational approach for bifurcation diagrams (period-doubling cascades, chaotic windows).

**Rich `ContinuousSystems.h` catalog:**  
Lorenz, Rössler, Van der Pol, Duffing, Hodgkin-Huxley, and others are fully implemented with literature-standard parameters. These serve as both pedagogical examples and validation targets. The Hodgkin-Huxley model (Nobel Prize-winning neuron model) in particular signals serious scientific ambition.

---

### 2.10 Functions Analyzer

`FunctionsAnalyzer.h` provides systematic function analysis:
- `SafeEvaluate()` — evaluates with exception catching, returns `SafeEvalResult{value, is_valid}`
- `FindZeros()` — bracket scan + Brent refinement
- `FindExtrema()` — sign changes in derivative
- `FindInflectionPoints()` — sign changes in second derivative
- Critical point classification: Maximum, Minimum, Saddle, Inflection, Asymptote, Discontinuity

The use of `SafeEvaluate` to handle functions with domain restrictions or divergences is a practical robustness measure — function analysis over broad domains often encounters discontinuities.

---

## 3. Weaknesses

### 3.1 Critical Issues

**W-C1: FFT limited to power-of-2 sizes**  
The Cooley-Tukey FFT requires `n = 2^k`. Signals of length 1000, 3000, or 10000 (common in instrumentation) cannot use the fast path and fall through to O(n²) DFT with no warning. In the worst case, this is a >10,000× performance regression for large n.  
**Fix:** Implement a mixed-radix FFT (Cooley-Tukey for arbitrary n via prime factorization). The Bluestein (chirp-z) algorithm as a fallback enables arbitrary n via power-of-2 transforms.

**W-C2: No gradient-free global optimization**  
`NelderMead` is a local method that converges to the nearest local minimum. `SimulatedAnnealing` and `GeneticAlgorithm` are referenced in AGENTS.md but not verified as implemented. True global optimization (multimodal functions, combinatorial-adjacent problems) requires methods like:
- Differential Evolution
- Basin hopping
- Random restart with multi-start NM

**Fix:** Verify SA/GA implementations; if absent, add Differential Evolution as the most widely applicable global derivative-free optimizer.

---

### 3.2 Significant Issues

**W-S1: `ODEStiffSolvers` — implementation status unclear**  
The ODE adaptive integrator provides `BulirschStoer_Stepper` and Dormand-Prince steppers for non-stiff systems. For stiff systems, the proper approach is implicit Runge-Kutta (Radau) or linearly-implicit methods. Whether the `ODEStiffSolvers.h` file provides complete, production-ready implementations or stubs needs verification. The presence of five DAE solvers suggests capability, but the bridge from "stiff ODE" to "DAE" requires the user to reformulate — which is non-trivial.  
**Fix:** Provide a clear stiff ODE solver (e.g., implicit trapezoidal / SDIRK) that does not require DAE formulation.

**W-S2: `BVPShootingMethod` — limited boundary condition types**  
Boundary value problems require flexible specification of boundary conditions. The shooting method handles separated BC (one set at left endpoint, one at right), but:
- Periodic BC: y(0) = y(T) (important for steady-state oscillations)
- Mixed/Robin BC: α·y + β·y' = γ at boundaries
- Integral constraints: ∫y dt = C

Without these, the BVP solver's usefulness is limited to simple test cases.  
**Fix:** Add periodic BC support and document the formulation for Robin BC.

**W-S3: `ChebyshevApproximation` — no adaptive degree selection**  
Chebyshev approximation requires specifying the polynomial degree in advance. For smooth functions with unknown smoothness properties, under-specifying yields large errors while over-specifying wastes computation. Adaptive degree selection (increase degree while ‖error‖ exceeds threshold) is the practically correct approach.  
**Fix:** Add `AdaptiveChebyshevApprox(f, a, b, tol)` that determines degree automatically by monitoring coefficients decay.

**W-S4: `Statistics` — no hypothesis testing beyond t-test**  
The statistics module supports mean, variance, moments, medians, and t-tests. Missing:
- **Chi-squared test** — goodness of fit, independence
- **Mann-Whitney U test** — non-parametric comparison
- **ANOVA** — multi-group comparison
- **Kolmogorov-Smirnov test** — distribution comparison
- **Correlation matrix** beyond simple `Correlation(x, y)`.

These are standard tools for any data analysis workflow.  
**Fix:** Add at minimum Chi-squared tests and one-way ANOVA.

**W-S5: `CompGeometry` — no 3D convex hull exposed in main API**  
`ConvexHull3D` exists but its exposure in the main include chain and documentation is thin relative to the 2D counterpart. For 3D mesh generation and physics simulations, 3D convex hull is critical.  
**Fix:** Ensure `ConvexHull3D` matches the `ConvexHull2D` API quality and document it prominently.

---

### 3.3 Minor Issues

**W-M1: `DynamicalSystemBase` — no default `integrate()` method**  
Every use of `DynamicalSystemBase<N,P>` requires the user to separately construct an `ODEAdaptiveIntegrator`, configure it, and pass the system. A convenience method:
```cpp
ODESystemSolution solve(const Vector<Real>& y0, Real t0, Real tend, Real eps = 1e-8)
```
on the base class would dramatically simplify the 95% case.

**W-M2: `DiscreteMaps` — missing Poincaré section integration**  
`DiscreteMaps.h` provides discrete maps (logistic, Hénon, standard). For continuous systems, Poincaré sections are obtained by running the ODE and sampling at discrete events — this is not automated. An `ODEPoincareSectionAnalyzer` that runs the ODE and records state when a condition is satisfied (e.g., `y[1] = 0, dy[1]/dt > 0`) would bridge continuous and discrete analysis.

**W-M3: `MatrixAlg` — matrix exponential method not documented**  
Matrix exponential (`exp(A)`) is critical for linear ODE solutions and control theory. Whether the implementation uses Padé approximation, Schur decomposition, or the scaling-and-squaring method affects accuracy for non-normal matrices. The method should be documented.

**W-M4: `GraphAlgorithms` — integration with Graph from base not described**  
`Graph.h` in the base layer defines the graph data structure. `GraphAlgorithms.h` presumably operates on it, but the relationship between `base/Graph.h`, `algorithms/GraphAlgorithms.h`, and `algorithms/GraphSpectral.h` is not clearly documented in the interface files. For users who want to use the graph tools, the entry point is unclear.

**W-M5: `FunctionsAnalyzer` — critical point bracketing may miss narrow features**  
The bracket-based scan for zeros and extrema uses a fixed step size `h`. For functions with narrow features (sharp peaks, rapid oscillations), `h` may skip over them entirely. The documentation should note sampling density as a user-tunable parameter and explain the tradeoff.

---

## 4. Algorithm Quality Benchmarking

### 4.1 Reference Comparison Table

| Category | MML | scipy/numpy | MATLAB | GSL |
|---|---|---|---|---|
| ODE Adaptive (non-stiff) | DP5, DP8, BS | RK45, DOP853 | ode45 | RK8PD |
| ODE Stiff | limited | Radau, BDF | ode15s | MSBDF |
| DAE | BDF2, BDF4, Radau, RODAS | ode15s (implicit only) | ode15i | none |
| Root Finding | Brent, Newton, Secant, Ridders | brentq, newton | fzero | gsl_root |
| Optimization 1D | Brent, GoldenSection | brent | fminbnd | gsl_min |
| Optimization ND | NM, Powell, CG, BFGS, LM | fmin, minimize | fminunc | multiroot |
| Eigenvalues | Jacobi, QR | eig (LAPACK) | eig | gsl_eigen |
| FFT | Cooley-Tukey radix-2 | FFTPACK/pocketfft | FFTW | gsl_fft |
| Statistics | basic | full scipy.stats | full stats toolbox | basic |
| Dynamical Systems | Lorenz++, Lyapunov, Bifurcation | none built-in | none built-in | none |

MML's DAE solver coverage and dynamical systems analysis tools are genuinely superior to scipy/GSL. The FFT limitation (power-of-2 only) is the clearest functional gap versus professional tools.

---

## 5. Systems Layer: Design Philosophy Assessment

The `mml/systems/` layer demonstrates a research-oriented mindset: Lyapunov exponents, bifurcation diagrams, and Hodgkin-Huxley neurons are topics from dynamical systems research, chaos theory, and computational neuroscience. This focus distinguishes MML as a library for **scientific computing research** rather than engineering applications. This is a strength (depth) and a limitation (breadth toward applied domains like optimization, control, and signal processing).

The `DynamicalSystemBase<N,P>` template design is forward-thinking: new systems are defined by subclassing and overriding `derivs()`. The base class handles all bookkeeping (dimension, parameter management, state naming) automatically. This is the correct object-oriented design for a "plug-in" system architecture.

---

## 6. Final Grade

| Category | Score (1–10) | Comments |
|---|---|---|
| **ODE solvers (non-stiff)** | 9.5 | FSAL, dense output, comprehensive stepper portfolio |
| **DAE solvers** | 9.0 | BDF2/4, Radau, RODAS — exceptional for single-header |
| **Root finding** | 9.0 | Brent, Newton+bracket, dual API, multi-root scanning |
| **Optimization** | 8.0 | Good 1D/ND coverage; global optimization status unclear |
| **Eigensolvers** | 8.5 | Jacobi correct, QR solid; LAPACK-style speed not achievable |
| **Statistics** | 7.0 | Core stats good; hypothesis testing limited |
| **Fourier** | 7.0 | DCT inclusion excellent; FFT power-of-2 only is a weak point |
| **Dynamical systems** | 9.5 | Lyapunov, bifurcation, fixed points — research-grade |
| **Comp. geometry** | 8.0 | KD-tree, Voronoi — solid; 3D consistently second-class |
| **Function analysis** | 8.5 | SafeEval, critical point detection well designed |

### Overall Algorithms & Systems Grade: **8.3 / 10**

The algorithms layer is where MML's ambition is most apparent. The DAE solver portfolio, Lyapunov analyzer, and bifurcation tools represent the kind of depth rarely seen outside specialized research libraries. The main weaknesses are the FFT restriction (power-of-2), unclear global optimization availability, and limited statistical hypothesis testing. The core ODE and root-finding capabilities are production-quality and match industry-standard tools.
