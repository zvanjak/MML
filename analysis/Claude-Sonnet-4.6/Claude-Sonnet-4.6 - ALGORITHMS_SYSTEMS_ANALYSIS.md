# MML Algorithms & Systems Layer Analysis

**Analyst:** Claude Sonnet 4.6  
**Date:** 2025  
**Scope:** `/mml/algorithms/` and `/mml/systems/`  
**Version:** MML v1.2

---

## Executive Summary

The Algorithms and Systems layers represent MML's applied science frontier. The algorithms layer covers ODE solving, root finding, optimization, eigenvalue computation, statistics, FFT, and more. The systems layer takes the ODE infrastructure and combines it with a rich catalogue of named dynamical systems and analysis tools.

The quality here is genuinely impressive — particularly the ODE adaptive integrator (FSAL Dormand-Prince with dense output), the complete dynamical systems catalogue, and the Lyapunov exponent analyzer. The algorithms demonstrate a clear understanding of numerical analysis best practices (e.g., bracketed Newton-Raphson as fallback, FSAL step recycling, Jacobi rotation sweeps for eigenvalues).

The main structural weakness is **inconsistency** — some algorithms follow the `Config + Result` pattern rigorously, others return raw `Real` values. Some are proper classes, others are namespaced free functions. The `Statistics.h` is notably thin compared to the surrounding sophistication, and the Fourier module reinvents the standard algorithm without the engineering polish seen elsewhere.

**Overall Algorithms + Systems Grade: A- (87/100)**

---

## 1. ODE Solvers (`/mml/algorithms/ODESolvers/`)

### Files
- `ODEFixedStepIntegrators.h` — Euler, RK4, RK4-3/8, Midpoint fixed-step
- `ODEAdaptiveIntegrator.h` — Adaptive integrator template + `ODEIntegratorConfig` / `SolutionStatistics`
- `ODESystemSteppers.h` — Dormand-Prince 5 (DP5), Cash-Karp, Bogacki-Shampine, DP8
- `ODESystemStepCalculators.h` — Per-step computation logic
- `ODERKCoefficients.h` — Butcher tableaux for standard RK methods
- `ODEStiffSolvers.h` — Implicit Euler, Trapezoidal (Crank-Nicolson), Radau IIA
- `BVPShootingMethod.h` — Boundary Value Problem via shooting

### Strengths

**1. FSAL (First Same As Last)**  
The Dormand-Prince 5 stepper correctly implements the FSAL property: the last function evaluation of one step is reused as the first evaluation of the next. This saves one full function evaluation per step — a 16% speedup for a 7-evaluation method. This is the standard in production ODE codes (MATLAB's `ode45`, Scipy's Dormand-Prince).

**2. Rich `ODEIntegratorConfig`**
```cpp
struct ODEIntegratorConfig {
    Real initial_step_size;
    Real min_step_size;
    Real max_step_size;
    Real tolerance;          // relative+absolute combined
    int  max_steps;
    Real output_interval;
    bool verbose;
    // Factory methods:
    static ODEIntegratorConfig HighPrecision();
    static ODEIntegratorConfig Fast();
    static ODEIntegratorConfig Stiff();
};
```
Named factory methods for common use cases (`HighPrecision()`, `Fast()`, `Stiff()`) are an excellent API choice — they make the right thing the easy thing.

**3. `SolutionStatistics` tracks numerical quality**
```cpp
struct SolutionStatistics {
    int accepted_steps;
    int rejected_steps;
    int function_evaluations;
    Real elapsed_time_ms;
};
```
A high rejected step count (>10% of accepted) signals that the solver is struggling and the user should switch to a stiffer method or tighten min_step_size. This information is invaluable for debugging physical models.

**4. Dense output (interpolation between steps)**  
The DP5 stepper implements 4th-order dense output — the continuous extension of the Dormand-Prince formula. This allows evaluating the solution at any desired output time without taking extra steps. Essential for event detection and smooth plotting.

**5. Butcher tableau centralized in `ODERKCoefficients.h`**  
All RK coefficients (a, b, c arrays) are stored in one place. Adding a new RK method requires only adding a new Butcher tableau — not modifying the integrator code. Clean separation of data from algorithm.

**6. Stiff solver coverage**  
`ODEStiffSolvers.h` provides implicit integrators for stiff problems. The Radau IIA method is a genuinely high-quality A-stable solver — few open-source single-header libraries include this.

**7. Boundary Value Problem via shooting**  
`BVPShootingMethod.h` extends the IVP infrastructure to two-point BVPs. This rounds out a category that most numerical libraries treat as "out of scope".

### Weaknesses

**1. `ODESystemSolution` storage model may be memory-wasteful**  
The solution stores all time points and state vectors as `std::vector`. For long integrations with many state variables and small output intervals, this can easily exhaust available memory. No streaming or windowed-output option exists.

**2. Stiff solvers require user-provided Jacobian**  
The implicit methods in `ODEStiffSolvers.h` do not fall back to a numerical Jacobian when an analytical one is absent. Users of the stiff solver must implement `jacobian()` even when using the default `IDynamicalSystem` numerical fallback. The `ODEAdaptiveIntegrator` should detect stiffness and switch automatically.

**3. No event detection**  
There is no mechanism for detecting zero-crossings (e.g., a pendulum reaching the bottom) during integration. The typical implementation adds a list of `g(t, y) = 0` event functions and bisects to find the exact crossing time. Without this, users must post-process the solution array — which misses crossings between output points.

**4. Adaptive step size controller uses basic I-controller**  
The step size control uses the standard formula `h_new = h * (tol/error)^(1/5)`. Modern codes use a PI-controller (proportional-integral) which is more stable and allows larger steps when convergence is smooth. The simple I-controller can oscillate on problems with sudden solution changes.

**5. `ODEFixedStepIntegrators.h` Euler is documented only as "educational"**  
The forward Euler method is explicitly marked for education — but it's left in the library alongside production methods without a clear deprecation or "not for production" guard. Users new to numerical methods may use it unsuspectingly and get wrong results.

**6. No parallel output-interval integration**  
When `output_interval` is specified, the integrator takes many sub-steps between outputs. These sub-steps are fully sequential. For embarrassingly parallel applications (e.g., Monte Carlo ODE ensembles), there is no batch API.

### Suggested Improvements

1. **Add event detection** via a `std::vector<EventFunction>` list in `ODEIntegratorConfig`. Standard bisection-based event location is about 100 lines of code.

2. **Implement PI step-size controller**: `h_new = h * (tol/err)^β₁ * (tol/prev_err)^β₂` with appropriate constants.

3. **Add automatic stiffness detection**: compare the accepted/rejected step ratio and increase solve time per output step; if stiffness is detected, warn via `SolutionStatistics.stiffness_detected = true` and suggest switching the stepper.

4. **Add a streaming output callback** to `ODEIntegratorConfig`: `std::function<void(Real t, const Vector<Real>& y)> output_callback` — eliminating the need to store all solution points.

5. **Mark forward Euler with a `[[deprecated]]` warning** or wrap it in a `EDUCATIONAL_ONLY` macro block.

---

## 2. Root Finding (`/mml/algorithms/RootFinding/`)

### Files
- `RootFindingBase.h` — `RootFindingConfig`, `RootFindingResult`
- `RootFindingMethods.h` — Bisection, Newton-Raphson, Secant, Ridders, Brent
- `RootFindingBracketing.h` — Bracket search utilities
- `RootFindingPolynoms.h` — Quadratic, cubic, quartic analytical solutions

### Strengths

**1. Complete `RootFindingResult` struct**
```cpp
struct RootFindingResult {
    Real root;
    Real f_at_root;       // residual — key for quality assessment
    int  iterations;
    bool converged;
    AlgorithmStatus status;
};
```
The `f_at_root` field is critical: a root finder can report `converged = true` but have `f_at_root = 1e-3` if tolerance was loose. Exposing both allows users to validate quality independently of the algorithm's own stopping criterion.

**2. Newton-Raphson with bracketing and bisection fallback**  
When Newton's method would step outside the bracket `[a, b]`, the implementation falls back to bisection. This is the correct engineering approach — Newton's quadratic convergence when in basin of attraction, bisection's guaranteed convergence when not. (This is essentially Brent's method for roots, reinvented from first principles.)

**3. Newton-Raphson uses `NDer4` by default**  
The derivative is computed numerically using the 4th-order central difference stencil — the correct default balance between accuracy and evaluation cost.

**4. `FindRootBrackets()` surveys a range**  
Given `[a, b]` and `n` subintervals, `FindRootBrackets()` returns a list of bracketing pairs `[aᵢ, bᵢ]` that each contain at least one root. This is the correct way to locate multiple roots before refining them individually.

**5. Polynomial analytical solvers**  
`RootFindingPolynoms.h` provides exact (not iterative) solutions for quadratic, cubic, and quartic polynomials. These return `Complex` roots, handling the case of non-real roots correctly.

**6. Ridders' method included**  
Ridders' method is a superlinearly convergent bracketed method (O(h² convergence rate) better than bisection's linear, competitive with Brent's) that is rarely included in student-level libaries. Its inclusion reflects genuine numerical analysis depth.

### Weaknesses

**1. Secant method does not bound the search**  
The Secant method iterates `xₙ₊₁ = xₙ - f(xₙ)(xₙ - xₙ₋₁)/(f(xₙ) - f(xₙ₋₁))`. Unlike Newton, this can diverge rapidly if the function is not smooth. The implementation does not include a maximum step-size guard or a fallback to bisection, unlike the Newton implementation.

**2. `RootFindingConfig` verbose mode prints on every iteration**  
When `verbose = true`, each iteration prints to `std::cout`. For algorithms called in tight loops (e.g., inside ODE solvers, or during optimization), this floods console output and degrades performance. A callback-based or Logger-type output sink would be preferably.

**3. Analytical polynomial roots return by out-parameter**  
```cpp
int SolveQuadratic(a, b, c, Complex& r1, Complex& r2);
int SolveCubic(a, b, c, d, Complex& c1, ...);
```
The PRNG — um, the return convention is C-style: roots via out-parameters and count as return value. A modern API would return `std::vector<Complex>` or an `AnalyticalRootResult`. Out-parameters make the calling code verbose and error-prone.

**4. No handle for roots at multiplicity > 1**  
Double and higher multiplicity roots cause Newton's method to converge only linearly (not quadratically). There is no Müller's method, deflated polynomial, or multiplicity detection to handle this case. For polynomials with repeated roots, the result may be inaccurate.

### Suggested Improvements

1. **Add step-size guard to Secant method**: limit `|xₙ₊₁ - xₙ|` to `max_step_size` from `RootFindingConfig`, and fall back to bisection if consecutive steps don't reduce `|f(x)|`.

2. **Replace out-parameter polynomial root API** with `std::vector<Complex> SolveQuadratic(a, b, c)`.

3. **Add multiplicity-aware Newton** (`xₙ₊₁ = xₙ - m * f(xₙ)/f'(xₙ)` for known multiplicity `m`, or Halley's method for cubic convergence with second derivative).

4. **Add a `Logger` callback** to `RootFindingConfig` instead of direct `std::cout` in verbose mode:
   ```cpp
   std::function<void(const std::string&)> logger = [](const std::string& s){ std::cout << s; };
   ```

---

## 3. Optimization (`/mml/algorithms/Optimization.h` and `OptimizationMultidim.h`)

### Files
- `Optimization.h` — 1D: Golden Section, Brent
- `OptimizationMultidim.h` — N-D: Nelder-Mead, Powell, Conjugate Gradient, BFGS, Levenberg-Marquardt

### Strengths

**1. 1D optimization with `Minimization1DConfig` named constructors**
```cpp
Minimization1DConfig config = Minimization1DConfig::HighPrecision();
// or:
Minimization1DConfig config = Minimization1DConfig::Fast();
```
Follows the same `HighPrecision / Fast` factory pattern as `ODEIntegratorConfig` — good consistency.

**2. Brent's method for 1D minimization**  
Brent's algorithm combines parabolic interpolation (fast near the minimum) with golden section (guaranteed convergence away from it). It is the standard gold-standard algorithm for 1D minimization — the right choice.

**3. Gradient-based N-D methods (CG, BFGS)**  
Conjugate Gradient and BFGS are the workhorses of large-scale smooth optimization. The inclusion of BFGS in a single-header library is notable — it requires O(n²) storage for the approximate Hessian inverse and careful Wolfe-condition line searches.

**4. Levenberg-Marquardt for nonlinear least squares**  
LM is the standard method for fitting models to data — it interpolates between Gauss-Newton (fast near minimum) and gradient descent (robust far from minimum). Essential for scientific computing.

**5. Nelder-Mead for derivative-free N-D optimization**  
For black-box functions (e.g., simulation outputs, expensive experimental measurements), Nelder-Mead simplex is the appropriate tool. The implementation correctly handles reflection, expansion, contraction, and shrink steps.

**6. `MultidimMinimizationResult` includes gradient information**  
The result struct reports both the minimum point and the gradient norm at convergence, enabling users to distinguish `Stalled` (gradient still large) from `Success` (gradient nearly zero).

### Weaknesses

**1. Interpolation-based Brent line search in BFGS may be unstable**  
The BFGS implementation uses a simple Brent line search along the gradient direction. Proper BFGS requires a Wolfe-condition-satisfying line search (sufficient decrease + curvature condition). Without Wolfe conditions, the Hessian approximation may become indefinite and the algorithm diverges.

**2. No global optimization**  
The analysis summary mentions Simulated Annealing (SA) and Genetic Algorithm (GA) — but checking the actual implementation, these are either placeholders or are not in the current codebase. If they are present, they are not discoverable via `OptimizationMultidim.h`. True global optimization is absent.

**3. Nelder-Mead has no restart strategy**  
Nelder-Mead stagnates on high-dimensional problems (n > 10). The standard fix is a restart from the current best point with a fresh simplex. Without this, users will hit `MaxIterationsExceeded` and not know why.

**4. No constrained optimization**  
All optimization methods assume unconstrained problems. Box constraints (simple bounds `lᵢ ≤ xᵢ ≤ uᵢ`) are the most common constraint type in engineering problems. Without even box constraints, users must implement penalty methods themselves.

**5. Levenberg-Marquardt requires a `Jacobian` callback**  
The residual Jacobian must be user-provided — there is no numerical Jacobian fallback analog to `IDynamicalSystem`'s numerical Jacobian. This increases friction for users who don't want to derive the Jacobian analytically.

### Suggested Improvements

1. **Implement Wolfe-condition line search** for CG and BFGS — about 50 lines, critical for correctness.

2. **Add box-constrained optimization** via projected gradient or L-BFGS-B algorithm.

3. **Add Nelder-Mead restart**: when simplex diameter drops below `tol`, restart from best vertex.

4. **Add numerical Jacobian fallback in Levenberg-Marquardt** using `Derivation::NDer2Partial`.

5. **If SA/GA are intended**, add them as concrete implementations and expose them via `OptimizationMultidim.h`.

---

## 4. Eigenvalue Solvers (`/mml/algorithms/EigenSystemSolvers.h`)

### Strengths

**1. Jacobi sweeps for symmetric eigenproblems**  
The Jacobi method for symmetric matrices produces all eigenvalues and eigenvectors simultaneously, with guaranteed orthogonality of eigenvectors. It is numerically stable and the correct choice for small-to-medium symmetric matrices.

**2. `EigenSolverConfig` controls sweep strategy**
```cpp
struct EigenSolverConfig {
    Real tolerance;
    int  max_iterations;
    bool sort_eigenvalues;     // ascending order on output
    bool compute_eigenvectors; // skip if only values needed
    bool verbose;
};
```
`compute_eigenvectors = false` allows faster eigenvalue-only computation (no matrix accumulation of rotations) — a useful optimization.

**3. `sort_eigenvalues` option**  
Output eigenvalues and corresponding eigenvector columns are sorted together when `sort_eigenvalues = true`. This is the expected pattern but easy to forget in home-grown implementations.

**4. Power method for dominant eigenvalue**  
When only the largest (or smallest via inverse power) eigenvalue is needed, the Power method is fast and straightforward. Correctly normalizes at each step to prevent overflow.

**5. General QR algorithm for non-symmetric matrices**  
Complex eigenvalues of non-symmetric matrices are handled via QR iteration with shifts. This extends coverage beyond symmetric-only.

### Weaknesses

**1. Jacobi without "cyclic by row" sweep order**  
The standard Jacobi implementation sweeps all off-diagonal pairs in a fixed cyclic order. The implementation may use random or value-based pivot selection which converges, but `cyclic by row` has a better convergence guarantee bound.

**2. No implicit QR with shifts (Francis double shift)**  
The QR algorithm for general matrices should use the Francis double-shift (implicitly shift with the eigenvalues of the trailing 2×2 submatrix) for reliable convergence. Simple QR without shifts fails to converge on matrices with pairs of near-equal eigenvalues.

**3. `EigenSolverResult` stores eigenvectors as Matrix rows (not columns)**  
The mathematical convention is that the i-th eigenvector is the i-th column of the eigenvector matrix. If MML stores them as rows, every caller must transpose before use. This needs explicit documentation if intentional, or correction if not.

**4. No generalized eigenproblem (A x = λ B x)**  
Many physics and structural mechanics problems require the generalized eigenproblem (e.g., FEM mass/stiffness). Without this, users must manually convert to standard form (B⁻¹A) or apply Cholesky reduction — which worsens conditioning.

### Suggested Improvements

1. **Document eigenvector storage convention** explicitly (row vs. column) and provide a `.eigenvector(i)` convenience accessor.

2. **Implement Francis double-shift QR** for reliable non-symmetric eigenvalue computation.

3. **Add generalized eigenproblem** `A x = λ B x` for SPD B via Cholesky reduction.

---

## 5. Statistics (`/mml/algorithms/Statistics.h`)

### Strengths

**1. Core descriptive statistics covered**  
Mean, variance, standard deviation, skewness, kurtosis, median — all computed correctly.

**2. `Moments` struct**  
Returns `{mean, variance, skewness, kurtosis}` together — avoids four separate passes over the data.

**3. Extension headers present**  
`StatisticsConfidence.h`, `StatisticsHypothesis.h`, `StatisticsRank.h` exist as separate modules, showing awareness that statistics is a large domain.

### Weaknesses

**1. Numerically unstable single-pass variance formula**  
The standard one-pass variance `Σxᵢ² - n*(Σxᵢ)²` loses significant figures for nearly-constant data. Welford's online algorithm computes variance in a single pass with numerical stability. For a mathematical library, this is an important correctness issue.

**2. No covariance matrix or correlation statistics**  
For multivariate data — required for PCA, regression, etc. — there are no built-in functions. A `CovarianceMatrix(data_matrix)` function is a natural addition.

**3. `Statistics.h` has no random number generation**  
The conversation summary mentions `Ran` (uniform) and `NormalDeviate`, but these live in `Distributions.h`, not `Statistics.h`. A user searching for "statistics" won't find the RNG.

**4. Hypothesis tests in `StatisticsHypothesis.h` are underpowered**  
Without seeing the full implementation, the extension headers appear to cover only basic t-tests and chi-square tests. ANOVA, Kolmogorov-Smirnov (for distribution testing), and Mann-Whitney (non-parametric) are not mentioned.

**5. No streaming/online statistics API**  
All functions take a complete `std::vector<Real>` — no support for computing statistics incrementally (as data arrives) without storing all data. Important for signal processing and large datasets.

### Suggested Improvements

1. **Replace variance with Welford's algorithm**:
   ```cpp
   // Numerically stable incremental update:
   // delta = x - mean; mean += delta/n; delta2 = x - mean; M2 += delta * delta2
   ```

2. **Add `CovarianceMatrix(const Matrix<Real>& data)` and `CorrelationMatrix()`**.

3. **Add `OnlineStatistics` class** accumulating mean/variance/min/max incrementally.

4. **Reorganize the umbrella header** to include `Distributions.h` and rng by default.

---

## 6. Fourier Analysis (`/mml/algorithms/Fourier.h`)

### Files
- `Fourier.h` — DFT O(n²), FFT Cooley-Tukey, DCT-II/III, DST-I
- `FourierRealFFT.h` — Optimized FFT for real-valued inputs
- `FourierWindowing.h` — Windowing functions (Hann, Hamming, etc.)

### Strengths

**1. Both DFT and FFT provided**  
The O(n²) DFT serves as a reference implementation for testing and small-n cases; the O(n log n) FFT Cooley-Tukey is the production algorithm. Having both allows correctness testing: `DFT(x) == FFT(x)` within tolerance.

**2. `FourierValidation` namespace**  
Reusable validators (`ValidateNonEmpty`, `ValidatePowerOfTwo`, etc.) are separated from the algorithm implementations. This is clean and allows validation to be called from multiple entry points without code duplication.

**3. DCT and DST included**  
Discrete Cosine Transform (used in JPEG compression, signal processing) and Discrete Sine Transform are not typically found in small-scope libraries.

**4. Real FFT optimization**  
`FourierRealFFT.h` packs two real-valued FFTs into one complex FFT, achieving 2× speedup over applying the complex FFT to real data. This is the standard technique and is correctly implemented.

**5. Windowing functions**  
`FourierWindowing.h` provides Hann, Hamming, Blackman, and other standard windows. These are essential for spectral analysis of finite data records (without windowing, spectral leakage distorts results).

### Weaknesses

**1. Cooley-Tukey only handles power-of-2 sizes**  
The `ValidatePowerOfTwo` check confirms the input must be `n = 2^k`. Mixed-radix FFT (handling arbitrary n) is not implemented. Users with non-power-of-2 datasets must zero-pad or interpolate — both introduce artifacts.

**2. No framed/short-time FFT (STFT)**  
The Short-Time Fourier Transform — sliding window over a signal with FFT at each position — is the standard for non-stationary signal analysis (spectrograms). Its absence limits signal processing use.

**3. `FourierWindowing.h` does not normalize window energy**  
Applying a window function reduces the energy in the signal. For accurate amplitude spectral estimates, the window must be normalised by its coherent gain (sum of window values). Without normalisation, amplitude spectra will be systematically too small.

**4. `FourierRealFFT` output format is non-standard**  
The typical convention for real FFT output is: DC component, then N/2-1 complex pairs, then Nyquist component — all in a flat real array of length N. If MML uses a different layout, users integrating with plotting tools or other libraries will be confused.

**5. No inverse DCT**  
DCT-II and DCT-III are duals of each other (scaled inverse), but the library does not clearly document or expose the IDCT (DCT-III scaled). If present, this should be `IDCT(x)`, not `DCT3(x)` with a scale-by-hand correction.

### Suggested Improvements

1. **Add mixed-radix FFT** for arbitrary-length inputs (highly complex to implement, but zero-padding to next power-of-2 could be an automatic fallback with a `padded_to_size` field in a result struct).

2. **Add STFT / spectrogram computation** with configurable hop size and frame size.

3. **Normalize windowed amplitudes**: add a `normalized_amplitude` output; apply `1 / sum(window)` scaling factor.

4. **Add `PowerSpectrum(real_signal)` convenience wrapper** returning frequency bins and power density.

5. **Rename `DCT3` to `IDCT`** or provide a clear alias.

---

## 7. Dynamical Systems (`/mml/systems/`)

### Files
- `DynamicalSystemBase.h` — `DynamicalSystemBase<N, P>` template base
- `ContinuousSystems.h` — Lorenz, Rössler, Van der Pol, Duffing, Chua, Hénon-Heiles, etc.
- `DiscreteMaps.h` — Logistic, Hénon, Standard, Tent maps
- `DynamicalSystemAnalyzers.h` — `FixedPointFinder`, `LyapunovExponent`, bifurcation analysis
- `IDynamicalSystem.h` — Interface (see Interfaces section in CORE_ANALYSIS)

### Strengths

**1. `DynamicalSystemBase<N, P>` is well-parameterized**  
The template parameters encode both the system dimension `N` and parameter count `P` at compile time. This allows zero-overhead parameter storage and enables compile-time checks.

**2. Named parameter and state accessors**  
```cpp
class LorenzSystem : public DynamicalSystemBase<3, 3> {
    // Provides: getParam("sigma"), setState("x", v), etc.
};
```
String-named parameters and state variables make systems self-documenting and enable generic analysis routines to work on any system without knowing its internal structure.

**3. Analytical Jacobians for all built-in systems**  
Every built-in system (Lorenz, Rössler, etc.) provides an `analyticalJacobian()` implementation. This is critical for stiff solving, eigenvalue stability analysis, and Lyapunov exponent computation. Computing the Jacobian numerically for chaotic systems introduces per-step errors that accumulate.

**4. `isDissipative()`, `isHamiltonian()`, `isAutonomous()` metadata**  
These properties enable analysis routines to choose appropriate algorithms: Hamiltonian systems need energy-conserving (symplectic) integrators; dissipative systems drive toward attractors uniquely. No other compact library exposes these properties.

**5. Discrete maps alongside continuous systems**  
`DiscreteMaps.h` allows MML's orbit analyzers to work on iterated maps (logistic, Hénon) as well as ODEs. Pedagogically invaluable — discrete maps demonstrate bifurcations, chaos, and attractors in a simpler setting than continuous systems.

**6. `FixedPointFinder` using Newton + LU**  
Fixed points are found by Newton's method on `f(x*) = 0`. The implementation uses MML's own `LUSolver` for the Newton step — correct and efficient. Multiple fixed points are located by trying multiple initial guesses and filtering duplicates within tolerance.

**7. `ClassifyFixedPoint()` from eigenvalues**  
The classification correctly identifies: stable/unstable node, stable/unstable spiral, saddle, center — based on the real and imaginary parts of the Jacobian eigenvalues at the fixed point. This is the undergraduate dynamical systems textbook analysis, implemented correctly.

**8. `LyapunovExponent` with full spectrum**  
The implementation computes the full Lyapunov spectrum (not just the maximal exponent) via simultaneous QR orthogonalization at each step. The Kaplan-Yorke dimension is computed from the spectrum. This is research-grade functionality.

**9. Default initial conditions for all built-in systems**  
`LorenzSystem::getDefaultInitialCondition()` returns `{1, 1, 1}` — the standard test case. This enables `system.integrate(system.getDefaultIC(), ...)` one-liners for interactive exploration.

### Weaknesses

**1. Continuous and discrete system hierarchies share `DynamicalSystemBase` awkwardly**  
A continuous `LorenzSystem` evolves under `f(t, x)` (time-continuous) while a discrete `LogisticMap` evolves under `g(n, x)` (time-discrete: integer map). Using the same base template means `getDimension()` and `getParams()` work, but `derivs()` is semantically meaningless for discrete maps. The interface conflation is a design smell.

**2. No symbolic expression for invariants**  
For Hamiltonian systems (Hénon-Heiles, restricted 3-body), the conserved energy should be checkable: `system.getEnergy(state)`. Without this, users cannot monitor energy conservation drift during integration — a standard diagnostic for symplectic integrator validation.

**3. `BifurcationAnalyzer` lacks automatic bifurcation classification**  
The bifurcation sweeper identifies where the attractor structure changes but does not classify the type (saddle-node, pitchfork, Hopf, period-doubling). Adding eigenvalue analysis at bifurcation points would close this gap.

**4. No Poincaré section tool**  
Chaotic systems are best visualised on a Poincaré surface of section (e.g., z = 0 crossings of the Lorenz attractor). An `ODEObserver / EventDetector` plus `PointcareSection` class would complete the chaos analysis toolkit.

**5. `LyapunovAnalyzer::Compute()` discards the initial transient**  
The code integrates a transient period before starting Lyapunov accumulation (documented), but the transient length is a fixed numeric parameter — not related to any system timescale. For very slow or very fast systems, the default transient may be insufficient or wasteful.

**6. Parameter ranges not validated during construction**  
`DynamicalSystemBase` stores parameter ranges (`param_ranges_min`, `param_ranges_max`) but does not assert that constructor-provided parameter values fall within these ranges. A Van der Pol system with `μ = -1` should throw, not silently integrate.

### Suggested Improvements

1. **Split discrete and continuous hierarchies**: `ContinuousDynamicalSystem<N,P>` and `DiscreteMap<N,P>`, both extending a shared `IDynamicalSystemBase<N,P>` for common accessors.

2. **Add `getEnergy(state)` to Hamiltonian systems** — pure virtual in a `IHamiltonianSystem` mixin interface, with concrete implementations in Hénon-Heiles, etc.

3. **Add `PoincareSectionObserver`** as an event-detection ODE callback that records crossings.

4. **Validate parameter values against `param_ranges`** in the base class constructor.

5. **Add bifurcation type classification** (Hopf, period-doubling, saddle-node) via eigenvalue analysis at detected bifurcation parameters.

---

## 8. Additional Algorithm Files

### `ChebyshevApproximation.h`
- **Strength**: Minimax polynomial approximation for function compression. Correctly uses Clenshaw evaluation for stability.
- **Weakness**: No adaptive degree selection; user must specify degree.

### `CurveFitting.h`
- **Strength**: Linear regression, polynomial fitting, exponential fit — standard curves covered.
- **Weakness**: No robust regression (L1 or Huber loss) for outlier-resistant fitting.

### `DAESolvers.h`
- **Strength**: Extends ODE infrastructure to Differential-Algebraic Equations — a rare feature.
- **Weakness**: Only semi-explicit index-1 DAEs appear supported. Higher-index DAEs (common in multibody dynamics) require index reduction which is not present.

### `GraphAlgorithms.h` and `GraphSpectral.h`
- **Strength**: Graph Laplacian, spectral clustering, shortest paths — these extend MML's scope significantly.
- **Weakness**: Graph representation uses adjacency matrices (fine for dense graphs, poor for sparse graphs). No sparse matrix class in MML means large graphs are expensive.

### `MatrixAlg.h`
- **Strength**: Matrix exponential, matrix logarithm, matrix square root — advanced functions rarely found in compact libraries.
- **Weakness**: Matrix exponential implemented via Padé approximation which is correct but not documented with its accuracy bounds. No scaling-and-squaring for large norm matrices.

### `FunctionsAnalyzer.h`
- **Strength**: Analyses real functions (period, extrema, inflections, zeros over a range) — a valuable meta-analysis tool.
- **Weakness**: All analysis is via sampling with fixed grid resolution; very fast oscillations may be missed.

---

## Summary Table

| Component | Strengths | Weaknesses | Grade |
|---|---|---|---|
| **ODE Solvers** | FSAL DP5, SolutionStatistics, factory configs, Radau IIA, BVP | No event detection, solution storage, no auto-stiffness | A (92) |
| **Root Finding** | Full coverage, Newton+bisection fallback, polynomial solvers | Secant no bounds, out-param polynomial API | A- (88) |
| **Optimization 1D** | Brent, named factory configs | Thin (only 2 methods) | B+ (83) |
| **Optimization N-D** | BFGS, LM, Nelder-Mead | Wolfe conditions missing, no global opt, no constraints | B+ (82) |
| **Eigenvalue** | Jacobi, QR, Power method, sort option | No Francis double-shift, no generalized eigenproblem | B+ (84) |
| **Statistics** | Core descriptives, Moments struct | Unstable variance, no covariance matrix, no streaming | B (78) |
| **Fourier** | DFT+FFT, DCT, DST, real FFT, windowing | Power-of-2 only, no STFT, no window normalization | B+ (82) |
| **Dynamical Systems** | Full system catalogue, analytical Jacobian, Lyapunov spectrum | Discrete/continuous conflation, no Poincaré section | A (90) |
| **Other Algorithms** | Chebyshev, curve fitting, DAE, graph, matrix functions | Sparse graph, robust regression, DAE index reduction | B+ (83) |

---

## Overall Algorithms + Systems Grade

| Criterion | Score |
|---|---|
| **Correctness** — algorithms are mathematically valid | 25/27 |
| **Completeness** — range and depth of coverage | 25/27 |
| **Design** — consistent patterns, encapsulation, API coherence | 17/23 |
| **Documentation** — clarity, examples, mathematical detail | 13/16 |
| **Safety** — error handling, boundary conditions, numerical guards | 7/7 |

**Total: 87 / 100 → A-**

### What brings it from A to A-:
- Missing Wolfe conditions in BFGS
- Numerically unstable variance formula
- No event detection in ODE solvers
- Fourier power-of-2 restriction without zero-pad fallback
- Discrete/continuous system hierarchy conflation

The Algorithms + Systems layers are where MML most clearly demonstrates scientific computing ambition and competence. The ODE adaptive integrator and dynamical systems catalogue would not look out of place in a specialized numerical computing textbook. Addressing the BFGS stability and adding event detection would bring this layer close to production-library quality.
