# Systems Package

**Dynamical Systems Analysis and the Unified Linear-Algebra Facade**

The `mml/systems/` layer provides two high-level facades built on MML's core primitives:

1. **`LinearSystem<Type>`** — one class that answers "solve and characterize my matrix": smart solver selection, decompositions, condition analysis, verification.
2. **The dynamical-systems framework** — fixed-point analysis, Lyapunov spectra, bifurcation diagrams, Poincaré sections, plus a library of canonical chaotic systems and discrete maps.

Models and analysis engines live in **`MML::Systems`**. Interfaces such as
`IDynamicalSystem` and `IDiscreteMap<N>`, and numerical primitives, live in **`MML`**.

## Features

### Dynamical Systems Analysis
- **Fixed Point Detection** — Newton–Raphson search (`FindFixedPoint`, `FindFixedPointsInBox`)
- **Stability Classification** — nodes, foci, saddles, centers (`FixedPointType`)
- **Lyapunov Exponents** — full spectrum via the Benettin QR method (`ComputeLyapunov`)
- **Bifurcation Diagrams** — parameter sweeps (`ComputeBifurcation`)
- **Poincaré Sections** — return maps and phase-space slicing (`ComputePoincareSection`)
- **One-call reports** — `Analyze()` produces a `DynamicalSystemReport` with a text summary

### Canonical Systems Included
- **Continuous** (`ContinuousSystems.h`): `LorenzSystem`, `RosslerSystem`, `VanDerPolSystem`, `DuffingSystem`, `ChuaCircuit`, `HenonHeilesSystem`, `DoublePendulumSystem` — all with analytic Jacobians where available
- **Discrete maps** (`DiscreteMaps.h`): `LogisticMap`, `HenonMap`, `StandardMap`, `TentMap`, plus `DiscreteMapLyapunov::Compute`

Built-in model constructors and parameter setters require finite coefficients, and rejected
bulk updates preserve all parameters. Double-pendulum masses and lengths must also be positive;
gravity may be finite and signed or zero. Suggested parameter ranges are advisory rather than
hard limits. Metadata indices are checked. See
[parameter domains and examples](DynamicalSystem.md#parameter-manipulation) and run
`MML_DocsApp --model-contracts` for the bounded executable examples.

### Linear Systems Interface
- **Unified Facade** — one class for all linear algebra on a given system
- **Smart Solver Selection** — `Solve()` picks triangular/Cholesky/QR/SVD/LU automatically
- **Matrix Decompositions** — `LUDecompose()`, `QRDecompose()`, `SVDDecompose()`, `CholeskyDecompose()`
- **Condition Analysis** — `ConditionNumber()`, `AssessStability()`, `ExpectedDigitsLost()`
- **Verification** — `Verify(x)` returns residuals and error estimates

## Quick Start

### Linear System Facade

The right-hand side is bound at construction; `Solve()` takes no arguments and auto-selects the best solver.

```cpp
#include <mml/systems/LinearSystem.h>

using namespace MML;
using namespace MML::Systems;

Matrix<Real> A(3, 3, { 4, 1, 2,
                       1, 5, 1,
                       2, 1, 6 });
Vector<Real> b({ 1, 2, 3 });

LinearSystem<Real> sys(A, b);

// Solve Ax = b (solver chosen automatically: SPD matrix -> Cholesky)
Vector<Real> x = sys.Solve();

// Verify the solution
auto check = sys.Verify(x);
std::cout << "Absolute residual: " << check.absoluteResidual << "\n";
std::cout << "Relative residual: " << check.relativeResidual << "\n";
std::cout << "Accurate: " << (check.isAccurate ? "yes" : "no") << "\n";

// Matrix properties (results cached per tolerance)
std::cout << "Symmetric:          " << sys.IsSymmetric() << "\n";
std::cout << "Positive definite:  " << sys.IsPositiveDefinite() << "\n";
std::cout << "Condition number:   " << sys.ConditionNumber() << "\n";
std::cout << "Determinant:        " << sys.Determinant() << "\n";

// Full analysis report
SystemAnalysis<Real> analysis = sys.Analyze();
std::cout << analysis.report;

// Decompositions (computed on demand, cached)
const auto& lu = sys.LUDecompose();    // lu.L, lu.U, lu.permutation, lu.determinant
const auto& qr = sys.QRDecompose();    // qr.Q, qr.R
const auto& svd = sys.SVDDecompose();  // svd.U, svd.V, svd.singularValues, svd.rank

// Eigenvalues (square matrices)
Vector<Complex> eigs = sys.Eigenvalues();

// Explicit solver choice when you know the structure
Vector<Real> xLU  = sys.SolveByLU();
Vector<Real> xChol = sys.SolveByCholesky();
Vector<Real> xLsq = sys.SolveLeastSquares();   // overdetermined systems (QR)

// Iterative solution for large sparse-ish systems
Vector<Real> xIter = sys.SolveIterative(IterativeMethod::Auto, 1e-10, 1000);
```

Multiple right-hand sides — factorize once, solve for every column:

```cpp
Matrix<Real> B(3, 2, { 1, 4,
                       2, 5,
                       3, 6 });
LinearSystem<Real> multi(A, B);
Matrix<Real> X = multi.SolveMultiple();   // A * X = B
```

There is also a one-line convenience free function:

```cpp
Vector<Real> x = Systems::SolveLinearSystem(A, b);
```

### Defining a Dynamical System

Derive from `DynamicalSystemBase<N, P>` (N = state dimension, P = number of parameters) and override `derivs`. Overriding `jacobian` is optional — a numerical Jacobian is used otherwise.

```cpp
#include <mml/systems/DynamicalSystem.h>   // umbrella header

using namespace MML;
using namespace MML::Systems;

// Most classic systems are already provided:
LorenzSystem lorenz;                    // sigma=10, rho=28, beta=8/3
VanDerPolSystem vdp(1.0);               // mu = 1.0

```

For a complete custom model, including validated constructor initialization and
a correctly sized analytical Jacobian, see the
[damped oscillator example](DynamicalSystem.md#idynamicalsystem-interface).
The bounded dynamical examples below run with `MML_DocsApp --systems-api`.
Use `DynamicalSystemAnalyzer<>` or `<Real>`: the engines use the selected `Real`
precision; the placeholder `Type` is not an independent mixed-precision switch.

### Fixed Point Analysis

```cpp
VanDerPolSystem vdp(1);
DynamicalSystemAnalyzer<> analyzer(vdp);
auto fixedPoint = analyzer.FindFixedPoint(Vector<Real>{0.1, 0.1});
if (fixedPoint.status != AlgorithmStatus::Success)
    throw std::runtime_error(fixedPoint.error_message);
std::cout << "Van der Pol linearized type: " << to_string(fixedPoint.type) << "\n";
auto points = analyzer.FindFixedPointsInBox(Vector<Real>{-2, -2}, Vector<Real>{2, 2}, 3);
std::cout << "Distinct successful box roots: " << points.size() << "\n";
```

These are non-exhaustive local searches for autonomous systems. A false
`isStable` does not prove instability; a linear `Center` does not prove nonlinear
stability or periodic motion. Use canonical `to_string`, not the deprecated `ToString`.

### Lyapunov Exponents

```cpp
LorenzSystem lorenz;
DynamicalSystemAnalyzer<> lorenzAnalyzer(lorenz);
Vector<Real> initial{1, 1, 1};
auto lyap = lorenzAnalyzer.ComputeLyapunov(initial, 20, 0.25, 0.01);
bool positiveGrowth = lorenzAnalyzer.IsChaotic(initial, 20);
Real growthTime = lorenzAnalyzer.GetLyapunovTime(initial, 20);
std::cout << "Finite-time spectrum: " << lyap.exponents << "\n";
std::cout << "Positive-growth heuristic: " << positiveGrowth << "\n";
std::cout << "Reciprocal detected growth rate: " << growthTime << "\n";
```

The three controls are total physical time, maximum QR interval and initial
adaptive IVP step, not a time step followed by an iteration count. Each quick
query performs a separate spectrum computation with its own defaults; reuse a
computed result when possible. A growth time is not a certified predictability
horizon. Check transient, duration and numerical-resolution sensitivity separately.

### Bifurcation Diagram

```cpp
LorenzSystem lorenz(10, 28, 8.0 / 3.0);
DynamicalSystemAnalyzer<> analyzer(lorenz);
auto bifurcation = analyzer.ComputeBifurcation(1, 24, 28, 3, Vector<Real>{1, 1, 1}, 2, 2, 4);
std::cout << "Sweep samples: " << bifurcation.parameterValues.size() << "\n";
```

This requests three parameter samples, records local maxima of component 2,
and uses two transient plus four recording time units per sample. The short
run illustrates the API, not an asymptotic bifurcation classification. State
and physical time continue between samples; the parameter is restored afterward.

For discrete maps use map iteration and `DiscreteMapLyapunov<N>` instead of
passing a map to the continuous `BifurcationAnalyzer::Sweep`:

```cpp
LogisticMap logistic(3.9);
auto orbit = logistic.orbit(Vector<Real>{0.3}, 32);
auto lyap = DiscreteMapLyapunov<1>::Compute(logistic, Vector<Real>{0.3}, 10000, 1000);
std::cout << "Logistic samples: " << orbit.size() << "\n";
std::cout << "Finite-iteration growth estimate: " << lyap.maxExponent << "\n";
```

The counts are measurement and discarded-transient iterations.
[Direct discrete parameter sampling](DynamicalSystem.md#discrete-parameter-samples)
is separate from a continuous-flow maxima sweep.

### Poincaré Section

```cpp
LorenzSystem lorenz;
DynamicalSystemAnalyzer<> analyzer(lorenz);
PoincareSection<Real> section(2, 27, +1);
auto crossings = analyzer.ComputePoincareSection(Vector<Real>{1, 1, 20}, section, 4, 0.01, 50, 1);
std::cout << "Bounded Poincare crossings: " << crossings.size() << "\n";
```

The search requests four full-state crossings of `z = 27` in the positive
direction, with a 50-time-unit search budget after a one-time-unit transient.
The simple API throws `ConvergenceError` on exhaustion;
`ComputePoincareSectionDetailed` instead retains partial crossings with
`MaxIterationsExceeded`.

### One-Call Full Analysis

```cpp
LorenzSystem lorenz;
DynamicalSystemAnalyzer<Real> analyzer(lorenz);

DynamicalSystemReport<Real> report = analyzer.Analyze(Vector<Real>({ 1, 1, 1 }), {}, 20);
std::cout << report.summary;   // local fixed points and qualified finite-time indicators
```

This bounded example runs with `MML_DocsApp --report-contracts`. The fixed-point
search is non-exhaustive; its stability counts describe local numerical
linearizations. Positive growth, a negative exponent sum, a near-zero sum, and
the Kaplan-Yorke estimate do not establish chaos, boundedness, Hamiltonian
structure, or a measured fractal dimension. The legacy heuristic flags retain
their existing numerical meanings. `Analyze` does not fill the trajectory or
Poincare collections; request those analyses separately.

## API Reference

### LinearSystem<Type>

Construction:

| Constructor | Description |
|-------------|-------------|
| `LinearSystem(A, b)` | System with a single RHS vector |
| `LinearSystem(A, B)` | System with multiple RHS columns |
| `LinearSystem(A)` | Analysis only, no RHS |

| Method | Description |
|--------|-------------|
| `Solve()` | Solve Ax = b, auto-selecting the solver |
| `SolveByLU() / SolveByCholesky() / SolveByQR() / SolveBySVD() / SolveByGaussJordan()` | Explicit solver choice |
| `SolveLeastSquares()` | Least-squares solution (QR) for overdetermined systems |
| `SolveIterative(method, tol, maxIter)` | Jacobi / Gauss–Seidel / SOR / auto |
| `SolveMultiple()` | Solve A·X = B for all columns (single factorization) |
| `Verify(x, tol)` | `VerificationResult`: residuals, backward/forward error |
| `Analyze()` | Full `SystemAnalysis` report |
| `IsSymmetric(tol) / IsPositiveDefinite(tol) / IsDiagonallyDominant()` | Structure queries (cached per tolerance) |
| `IsUpperTriangular(tol) / IsLowerTriangular(tol) / IsDiagonal(tol)` | Shape queries |
| `Determinant() / Rank() / Nullity()` | Basic invariants |
| `ConditionNumber() / ConditionNumber1() / ConditionNumberInfinity()` | Condition estimates |
| `AssessStability() / ExpectedDigitsLost()` | Numerical health |
| `LUDecompose() / QRDecompose() / SVDDecompose() / CholeskyDecompose()` | Cached decompositions |
| `NullSpace() / ColumnSpace() / RowSpace() / LeftNullSpace()` | Fundamental subspaces (SVD-based) |
| `Inverse() / PseudoInverse()` | Matrix inverses |
| `Eigensystem() / Eigenvalues() / SymmetricEigenvalues()` | Eigen analysis |
| `SpectralRadius() / HasComplexEigenvalues()` | Spectral queries |

### DynamicalSystemAnalyzer<Type>  (Type is a placeholder: Real only)

Constructed from a reference to any `IDynamicalSystem` implementation:

| Method | Description |
|--------|-------------|
| `FindFixedPoint(guess, tol)` | Newton search from one guess |
| `FindFixedPoints(guesses)` | Newton search from many guesses (deduplicated) |
| `FindFixedPointsInBox(min, max, gridPerDim)` | Grid search over a box |
| `ComputeLyapunov(x0, totalTime)` | Full Lyapunov spectrum (Benettin QR) |
| `ComputeMaxLyapunov(x0, totalTime)` | Largest entry from a full-spectrum computation |
| `ComputeBifurcation(paramIdx, pMin, pMax, nSteps, x0)` | Parameter sweep |
| `ComputePoincareSection(x0, section, nCrossings, h, maxSearchTime, transientTime)` | Section crossings (bounded search; throws `ConvergenceError` if the budget is exhausted) |
| `IntegrateTrajectory(x0, totalTime)` | Raw trajectory points |
| `IsChaotic(x0) / IsDissipative(x0)` | Legacy finite-time positive-growth / negative-sum heuristics, not global classifications |
| `GetLyapunovTime(x0) / GetFractalDimension(x0)` | Reciprocal detected growth rate / spectrum-based Kaplan-Yorke estimate |
| `ScanChaosVsParameter(paramIdx, pMin, pMax, nSteps, x0)` | (param, λ₁) pairs |
| `Analyze(x0)` | Local fixed-point searches and finite-time Lyapunov indicators; read `report.summary` for qualified text |

Lower-level static engines (used by the facade, callable directly):
`FixedPointFinder::Find/FindMultiple`, `LyapunovAnalyzer::Compute`,
`BifurcationAnalyzer::Sweep`, `PhaseSpaceAnalyzer::ComputePoincareSection/IntegrateTrajectory`.

### IDynamicalSystem / DynamicalSystemBase<N, P>

| Method | Description |
|--------|-------------|
| `getDim()` | State-space dimension |
| `derivs(t, y, dydt)` | Compute dy/dt = f(t, y) — **must override** |
| `jacobian(t, y, J)` | Jacobian matrix (numerical default) |
| `hasAnalyticalJacobian()` | Capability declaration; overrides must advertise it accurately |
| `getNumParam() / getParam(i) / setParam(i, v)` | Parameter access |
| `getStateName(i) / getParamName(i)` | Naming for reports |
| `isAutonomous() / isHamiltonian() / isDissipative()` | System character flags |
| `getDefaultInitialCondition()` | Model default; inherited zero may be an equilibrium |
| `getNumInvariants() / computeInvariant(i, x)` | Conserved quantities (e.g. energy) |

### Fixed Point Types

| Type | Eigenvalue Pattern | Behavior |
|------|-------------------|----------|
| `StableNode` | All λ < 0 (real) | Trajectories converge |
| `UnstableNode` | All λ > 0 (real) | Trajectories diverge |
| `Saddle` | Mixed signs | Stable/unstable manifolds |
| `StableFocus` | Re(λ) < 0 (complex) | Spiral inward |
| `UnstableFocus` | Re(λ) > 0 (complex) | Spiral outward |
| `Center` | Nonzero purely imaginary spectrum within the tolerance band | Linearized label only; no nonlinear stability or periodic-orbit proof |

### Result Structures (DynamicalSystemTypes.h)

These library types are already declared by the umbrella header; do not redeclare
partial lookalikes in client code.

| Type | Important fields and interpretation |
|------|-------------------------------------|
| `FixedPoint<Real>` | `location`, `jacobian`, `eigenvalues`, `type`, `isStable`, `convergenceResidual`, `iterations`, `status`, `error_message`, `function_evaluations`; check status before treating a candidate as a classified equilibrium |
| `LyapunovResult<Real>` | `exponents`, `maxExponent`, `sum`, `kaplanYorkeDimension`, `isChaotic`, `numOrthonormalizations`, `totalTime`, `function_evaluations`, `qr_refinements`; finite-time estimates and work diagnostics, not chaos or dimension certificates |
| `BifurcationDiagram<Real>` | `parameterName`, `parameterValues`, `attractorValues`; legacy `numTransientSteps` and `numRecordedPoints` are rounded phase-duration/initial-step estimates, not actual solver work or stored counts |
| `PoincareSection<Real>` | `variable`, `value`, `direction`; direction is -1, 0, or +1; crossing APIs return full states |

See [Detailed evaluation contracts](DynamicalSystem.md#detailed-evaluation-contracts)
for mandatory safety checks, status conversion, actual direct RHS counts, partial
results, and validated zero-work requests.

## File Structure

```
mml/systems/
├── LinearSystem.h              # Unified linear-algebra facade
├── DynamicalSystem.h           # Umbrella header (includes all below)
├── DynamicalSystem/
│   ├── DynamicalSystemTypes.h  # FixedPoint, LyapunovResult, BifurcationDiagram, PoincareSection
│   ├── DynamicalSystemBase.h   # DynamicalSystemBase<N, P> base class
│   ├── DynamicalAnalysisCommon.h # Shared integration and detailed-API support
│   ├── DynamicalSystemValidation.h # Shared validation and explicit parameter restoration
│   ├── LyapunovQR.h             # Scaled tangent-frame orthonormalization
│   ├── FixedPointAnalysis.h    # FixedPointFinder
│   ├── LyapunovAnalysis.h      # LyapunovAnalyzer
│   ├── BifurcationAnalysis.h   # BifurcationAnalyzer
│   └── PhaseSpaceAnalysis.h    # PhaseSpaceAnalyzer
├── DynamicalSystemAnalyzer.h   # DynamicalSystemAnalyzer facade + DynamicalSystemReport
├── ContinuousSystems.h         # Lorenz, Rössler, Van der Pol, Duffing, Chua, Hénon-Heiles, double pendulum
└── DiscreteMaps.h              # Logistic, Hénon, Standard, Tent maps + DiscreteMapLyapunov
```

Tests live in `tests/systems/`.

## Mathematical Background

### Lyapunov Exponents

The Lyapunov exponent measures the rate of separation of infinitesimally close trajectories:

$$\lambda = \lim_{t \to \infty} \frac{1}{t} \ln \frac{\|\delta x(t)\|}{\|\delta x(0)\|}$$

Computed via repeated QR re-orthonormalization of the variational equation solution (Benettin's method).

### Fixed Point Stability

At equilibrium $x^*$ where $f(x^*) = 0$, linearize:

$$\dot{y} = Df(x^*) \cdot y$$

For sufficiently smooth autonomous systems, strictly negative real parts imply
local asymptotic stability and a positive real part implies instability.
Neutral eigenvalues require further analysis; a `Center` label alone proves
neither nonlinear stability nor periodicity.

### Kaplan-Yorke Dimension

$$D_{KY} = j + \frac{\sum_{i=1}^{j} \lambda_i}{|\lambda_{j+1}|}$$

where $j$ is largest integer such that $\sum_{i=1}^{j} \lambda_i \geq 0$.
The interpolation formula applies when a next exponent exists and its denominator
is nonzero. Using a finite-time spectrum yields a heuristic dimension estimate;
the analysis does not verify the assumptions identifying it with an invariant
set's dimension.

### Condition Number

$$\kappa(A) = \|A\| \cdot \|A^{-1}\| = \frac{\sigma_{\max}}{\sigma_{\min}}$$

Measures sensitivity of linear system solution to perturbations.

## Included Example Systems

### Classic Chaotic Systems

| Class | Dimension | Behavior |
|-------|-----------|----------|
| `LorenzSystem` | 3D | Strange attractor, chaos |
| `RosslerSystem` | 3D | Simpler chaotic attractor |
| `ChuaCircuit` | 3D | Double-scroll attractor |
| `HenonMap` | 2D (map) | Fractal structure |
| `LogisticMap` | 1D (map) | Period-doubling to chaos |
| `TentMap` | 1D (map) | Piecewise-linear chaos |
| `StandardMap` | 2D (map) | Area-preserving, KAM tori |

### Classic Oscillators and Conservative Systems

| Class | Dimension | Behavior |
|-------|-----------|----------|
| `VanDerPolSystem` | 2D | Limit cycle |
| `DuffingSystem` | 3D (forced) | Nonlinear oscillator |
| `HenonHeilesSystem` | 4D | Hamiltonian, mixed phase space |
| `DoublePendulumSystem` | 4D | Conservative chaos |

## References

- Strogatz, S.H. (2015). "Nonlinear Dynamics and Chaos" (2nd ed.)
- Ott, E. (2002). "Chaos in Dynamical Systems" (2nd ed.)
- Parker, T.S., and Chua, L.O. (1989). "Practical Numerical Algorithms for Chaotic Systems"
- Numerical Recipes in C++, 3rd Edition, Chapter 17: Integration of ODEs

## See Also

- [mml/algorithms/](../../mml/algorithms/) — ODE solvers used for trajectory work
- [mml/core/LinAlgEqSolvers/](../../mml/core/LinAlgEqSolvers/) — the solver engines behind `LinearSystem`
- [tests/systems/](../../tests/systems/) — usage examples in test form
