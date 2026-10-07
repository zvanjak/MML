# Systems Package

**Dynamical Systems Analysis and the Unified Linear-Algebra Facade**

The `mml/systems/` layer provides two high-level facades built on MML's core primitives:

1. **`LinearSystem<Type>`** — one class that answers "solve and characterize my matrix": smart solver selection, decompositions, condition analysis, verification.
2. **The dynamical-systems framework** — fixed-point analysis, Lyapunov spectra, bifurcation diagrams, Poincaré sections, plus a library of canonical chaotic systems and discrete maps.

Everything in this package lives in the **`MML::Systems`** namespace.

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

// Or define your own:
class MySystem : public DynamicalSystemBase<2, 1> {
public:
    MySystem(Real mu = 1.0) {
        _params[0] = mu;
        _stateNames = { "x", "y" };
        _paramNames = { "mu" };
    }

    void derivs(Real /*t*/, const Vector<Real>& y, Vector<Real>& dydt) const override {
        Real mu = _params[0];
        dydt[0] = y[1];
        dydt[1] = mu * (1 - y[0] * y[0]) * y[1] - y[0];
    }
};
```

### Fixed Point Analysis

```cpp
DynamicalSystemAnalyzer<Real> analyzer(vdp);

// Single fixed point from an initial guess
FixedPoint<Real> fp = analyzer.FindFixedPoint(Vector<Real>({ 0.1, 0.1 }));
std::cout << "Location: " << fp.location << "\n";
std::cout << "Type:     " << ToString(fp.type) << "\n";
std::cout << "Stable:   " << fp.isStable << "\n";

// Search a whole box on a grid of starting guesses
auto points = analyzer.FindFixedPointsInBox(Vector<Real>({ -2, -2 }),
                                            Vector<Real>({  2,  2 }), 5);
for (const auto& p : points)
    std::cout << p.location << "  " << ToString(p.type) << "\n";
```

### Lyapunov Exponents

```cpp
LorenzSystem lorenz;
DynamicalSystemAnalyzer<Real> analyzer(lorenz);

Vector<Real> x0({ 1.0, 1.0, 1.0 });
LyapunovResult<Real> lyap = analyzer.ComputeLyapunov(x0, 1000.0);

std::cout << "Exponents:      " << lyap.exponents << "\n";
std::cout << "Max (lambda_1): " << lyap.maxExponent << "\n";
std::cout << "Sum:            " << lyap.sum << "\n";           // < 0 for dissipative
std::cout << "K-Y dimension:  " << lyap.kaplanYorkeDimension << "\n";
std::cout << "Chaotic:        " << (lyap.isChaotic ? "yes" : "no") << "\n";

// Quick queries
bool chaotic = analyzer.IsChaotic(x0);
Real tPredict = analyzer.GetLyapunovTime(x0);   // 1 / lambda_max
```

### Bifurcation Diagram

```cpp
LorenzSystem lorenz;
DynamicalSystemAnalyzer<Real> analyzer(lorenz);

// Sweep parameter 1 (rho) over [20, 30] in 100 steps
BifurcationDiagram<Real> diagram =
    analyzer.ComputeBifurcation(1, 20.0, 30.0, 100, Vector<Real>({ 0.1, 0.1, 0.1 }));

// diagram.parameterValues[i] and diagram.attractorValues[i] (local maxima of a
// state component after transients) can be plotted directly.
```

For discrete maps use the map classes and `DiscreteMapLyapunov`:

```cpp
LogisticMap logistic(4.0);
auto orbit = logistic.orbit(Vector<Real>({ 0.3 }), 1000);
auto lyapMap = DiscreteMapLyapunov::Compute(logistic, Vector<Real>({ 0.3 }));
// LogisticMap at r=4 has analytical Lyapunov exponent ln(2)
```

### Poincaré Section

```cpp
LorenzSystem lorenz;
DynamicalSystemAnalyzer<Real> analyzer(lorenz);

Vector<Real> x0({ 1.0, 1.0, 20.0 });

// Section: state variable 2 (z) crossing z = 27 in the positive direction
PoincareSection<Real> section(2, 27.0, +1);

std::vector<Vector<Real>> crossings =
    analyzer.ComputePoincareSection(x0, section, 1000);

// Each element is the full state at a crossing - plot (x, y) to see the
// characteristic Lorenz attractor cross-section.
```

### One-Call Full Analysis

```cpp
LorenzSystem lorenz;
DynamicalSystemAnalyzer<Real> analyzer(lorenz);

DynamicalSystemReport<Real> report = analyzer.Analyze(Vector<Real>({ 1, 1, 1 }));
std::cout << report.summary;   // fixed points, Lyapunov spectrum, chaos verdict
```

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

### DynamicalSystemAnalyzer<Type>

Constructed from a reference to any `IDynamicalSystem` implementation:

| Method | Description |
|--------|-------------|
| `FindFixedPoint(guess, tol)` | Newton search from one guess |
| `FindFixedPoints(guesses)` | Newton search from many guesses (deduplicated) |
| `FindFixedPointsInBox(min, max, gridPerDim)` | Grid search over a box |
| `ComputeLyapunov(x0, totalTime)` | Full Lyapunov spectrum (Benettin QR) |
| `ComputeMaxLyapunov(x0, totalTime)` | Largest exponent only |
| `ComputeBifurcation(paramIdx, pMin, pMax, nSteps, x0)` | Parameter sweep |
| `ComputePoincareSection(x0, section, nCrossings)` | Section crossings |
| `IntegrateTrajectory(x0, totalTime)` | Raw trajectory points |
| `IsChaotic(x0) / IsDissipative(x0)` | Boolean chaos/dissipation tests |
| `GetLyapunovTime(x0) / GetFractalDimension(x0)` | Derived scalar measures |
| `ScanChaosVsParameter(paramIdx, pMin, pMax, nSteps, x0)` | (param, λ₁) pairs |
| `Analyze(x0)` | Everything above in one `DynamicalSystemReport` |
| `GenerateSummary(report)` | Human-readable text summary |

Lower-level static engines (used by the facade, callable directly):
`FixedPointFinder::Find/FindMultiple`, `LyapunovAnalyzer::Compute`,
`BifurcationAnalyzer::Sweep`, `PhaseSpaceAnalyzer::ComputePoincareSection/IntegrateTrajectory`.

### IDynamicalSystem / DynamicalSystemBase<N, P>

| Method | Description |
|--------|-------------|
| `getDim()` | State-space dimension |
| `derivs(t, y, dydt)` | Compute dy/dt = f(t, y) — **must override** |
| `jacobian(t, y, J)` | Jacobian matrix (numerical default) |
| `hasAnalyticalJacobian()` | True if `jacobian` is exact |
| `getNumParam() / getParam(i) / setParam(i, v)` | Parameter access |
| `getStateName(i) / getParamName(i)` | Naming for reports |
| `isAutonomous() / isHamiltonian() / isDissipative()` | System character flags |
| `getDefaultInitialCondition()` | Sensible starting state |
| `getNumInvariants() / computeInvariant(i, x)` | Conserved quantities (e.g. energy) |

### Fixed Point Types

| Type | Eigenvalue Pattern | Behavior |
|------|-------------------|----------|
| `StableNode` | All λ < 0 (real) | Trajectories converge |
| `UnstableNode` | All λ > 0 (real) | Trajectories diverge |
| `Saddle` | Mixed signs | Stable/unstable manifolds |
| `StableFocus` | Re(λ) < 0 (complex) | Spiral inward |
| `UnstableFocus` | Re(λ) > 0 (complex) | Spiral outward |
| `Center` | Re(λ) = 0 (imaginary) | Periodic orbits |

### Result Structures (DynamicalSystemTypes.h)

```cpp
template<typename Type = Real>
struct FixedPoint {
    Vector<Type> location;                         // Position in state space
    std::vector<std::complex<Type>> eigenvalues;   // Jacobian eigenvalues
    Matrix<Type> jacobian;                         // Jacobian at fixed point
    FixedPointType type;                           // Stability classification
    bool isStable;                                 // Overall stability
    Type convergenceResidual;                      // ||f(x*)|| at convergence
    int iterations;                                // Newton iterations used
};

template<typename Type = Real>
struct LyapunovResult {
    Vector<Type> exponents;         // λ₁ ≥ λ₂ ≥ ... ≥ λₙ
    Type maxExponent;               // λ₁
    Type sum;                       // Σλᵢ
    Type kaplanYorkeDimension;      // Fractal dimension estimate
    bool isChaotic;                 // True if λ₁ > 0
    int numOrthonormalizations;     // QR steps performed
    Type totalTime;                 // Integration time
};

template<typename Type = Real>
struct BifurcationDiagram {
    std::string parameterName;
    std::vector<Type> parameterValues;
    std::vector<std::vector<Type>> attractorValues;
    int numTransientSteps;
    int numRecordedPoints;
};

template<typename Type = Real>
struct PoincareSection {
    int variable;    // Which state variable defines the section
    Type value;      // Section at x[variable] = value
    int direction;   // +1 positive crossing, -1 negative, 0 both

    PoincareSection(int var = 0, Type val = 0, int dir = 0);
};
```

## File Structure

```
mml/systems/
├── README.md                   # This file
├── LinearSystem.h              # Unified linear-algebra facade
├── DynamicalSystem.h           # Umbrella header (includes all below)
├── DynamicalSystem/
│   ├── DynamicalSystemTypes.h  # FixedPoint, LyapunovResult, BifurcationDiagram, PoincareSection
│   ├── DynamicalSystemBase.h   # DynamicalSystemBase<N, P> base class
│   ├── DynamicalAnalysisCommon.h # Shared integration and detailed-API support
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

Stability determined by eigenvalues of Jacobian $Df(x^*)$.

### Kaplan-Yorke Dimension

$$D_{KY} = j + \frac{\sum_{i=1}^{j} \lambda_i}{|\lambda_{j+1}|}$$

where $j$ is largest integer such that $\sum_{i=1}^{j} \lambda_i \geq 0$.

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

- [mml/algorithms/](../algorithms/) — ODE solvers used for high-accuracy trajectory work
- [mml/core/LinAlgEqSolvers/](../core/LinAlgEqSolvers/) — the solver engines behind `LinearSystem`
- [tests/systems/](../../tests/systems/) — usage examples in test form
