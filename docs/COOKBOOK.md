# MML Cookbook

**Task-oriented recipes** — each one takes a concrete problem and shows the complete path
from setup to solution. For the library-wide ideas the recipes build on (interfaces,
Config + Result, contracts), read [Fundamentals.md](Fundamentals.md) first.

**Recipe format** (fixed):
**1. Setup** — state the problem, construct the inputs ·
**2. Solve** — call MML, interpret the result (convergence, diagnostics).
Per the docs↔demos contract (AGENTS.md), **every recipe is backed by a compilable,
runnable demo** in [src/docs_demos/docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp)
(one file for all recipes, target `MML_DocsApp`) — a recipe without its demo is not done.

## High-Level Table of Contents

1. [Solving linear systems](#recipe-1-solving-linear-systems)
2. [Root finding](#recipe-2-root-finding)
3. [Calculating integrals](#recipe-3-calculating-integrals)
4. [Calculating derivatives](#recipe-4-calculating-derivatives)
5. [Optimization recipes](#recipe-5-optimization-recipes)
6. [Solving differential equations](#recipe-6-solving-differential-equations)
7. [Eigenvalues and eigenvectors](#recipe-7-eigenvalues-and-eigenvectors)
8. [Matrix decompositions](#recipe-8-matrix-decompositions)
9. [Matrix properties and diagnostics](#recipe-9-matrix-properties-and-diagnostics)
10. [Statistics](#recipe-10-statistics)
11. [3D geometry: planes, lines, and bodies](#recipe-11-3d-geometry-planes-lines-and-bodies)
12. [3D geometry: algorithms](#recipe-12-3d-geometry-algorithms)
13. [Interpolating functions](#recipe-13-interpolating-functions)
14. [Curve fitting](#recipe-14-curve-fitting-including-fitting-to-points)
15. [Curves and surfaces](#recipe-15-curves-and-surfaces)
16. [Polynoms: polynomial models and approximations](#recipe-16-polynoms-polynomial-models-and-approximations)
17. [Quaternions](#recipe-17-quaternions)
18. [Analyzing functions](#recipe-18-analyzing-functions)
19. [Coordinate transformations](#recipe-19-coordinate-transformations)
20. [Fields](#recipe-20-fields)
21. [Tensors](#recipe-21-tensors)
22. [Differential geometry and geodesics](#recipe-22-differential-geometry-and-geodesics)
23. [Fourier algorithms, including FFT](#recipe-23-fourier-algorithms)
24. [Graph algorithms](#recipe-24-graph-algorithms)
25. [DAE solvers](#recipe-25-dae-solvers)
26. [Algebra: finite structures, symmetry, and exact arithmetic](#recipe-26-algebra-finite-structures-symmetry-and-exact-arithmetic)
27. [Complex analysis](#recipe-27-complex-analysis)
28. Function spaces *(planned)*
29. Vector spaces *(planned)*
30. [Combinatorics and number theory](#recipe-30-combinatorics-and-number-theory)
31. Using special functions *(planned)*
32. Persistence recipes *(planned)*
33. Serialization recipes *(planned)*
34. Visualization recipes *(planned)*


---

## Detailed Table of Contents

### Basic Recipes

1. ✅ **[Solving linear systems](#recipe-1-solving-linear-systems)**
   - [1.1 Real 3×3 system — LU direct, automatic selection, `LinearSystem` facade](#11-real-33-system)
   - [1.2 Complex 3×3 system](#12-complex-33-system)
2. ✅ **[Root finding](#recipe-2-root-finding)**
   - [2.1 Classic: bracket the root, then solve (Kepler's equation)](#21-classic-bracket-the-root-then-solve)
   - [2.2 All roots in an interval (beam frequencies)](#22-all-roots-in-an-interval)
   - [2.3 Polynomial roots (Wilkinson-5)](#23-polynomial-roots)
   - [2.4 Roots of a complex function (eᶻ = z)](#24-roots-of-a-complex-function)
   - [2.5 Roots in multiple dimensions (nonlinear system)](#25-roots-in-multiple-dimensions)
3. ✅ **[Calculating integrals](#recipe-3-calculating-integrals)**
   - [3.1 Real function — trapezoid, Romberg, Gauss-Kronrod](#31-real-function--trapezoid-romberg-gauss-kronrod)
   - [3.2 2D integration — simple and adaptive](#32-2d-integration--simple-and-adaptive)
   - [3.3 3D integration — simple and adaptive](#33-3d-integration--simple-and-adaptive)
   - [3.4 Improper integrals — infinite ranges, endpoint singularity](#34-improper-integrals)
   - [3.5 Path integration — work done in a force field](#35-path-integration--work-done-in-a-force-field)
   - [3.6 Surface integration — flux through a closed surface](#36-surface-integration--flux-through-a-closed-surface)
   - [3.7 Monte Carlo integration](#37-monte-carlo-integration)
4. ✅ **[Calculating derivatives](#recipe-4-calculating-derivatives)**
   - [4.1 Real function — first, second, and third derivatives](#41-real-function--first-second-and-third-derivatives)
   - [4.2 Scalar function — one partial and all partials](#42-scalar-function--one-partial-and-all-partials)
   - [4.3 Vector function — partial, partial-by-all, all-by-all](#43-vector-function--partial-partial-by-all-all-by-all)
   - [4.4 Parametric curve — first, second, and third derivatives](#44-parametric-curve--first-second-and-third-derivatives)
   - [4.5 Parametric surface — tangent plane, normal, and area scale](#45-parametric-surface--tangent-plane-normal-and-area-scale)
   - [4.6 Jacobians and Hessians](#46-jacobians-and-hessians)
   - [4.7 High-precision complex-step differentiation](#47-high-precision-complex-step-differentiation)
5. ✅ **[Optimization recipes](#recipe-5-optimization-recipes)**
    - [5.1 Real function minimum finding — golden section and Brent](#51-real-function--golden-section-and-brent)
    - [5.2 Multidimensional function optimization — Powell, Nelder-Mead, and BFGS](#52-multidimensional-function--powell-nelder-mead-and-bfgs)
    - [5.3 Linear programming example](#53-linear-programming-example)
6. ✅ **[Solving differential equations](#recipe-6-solving-differential-equations)**
   - [6.1 Fixed-step RK4 on a simple problem](#61-fixed-step-rk4-on-a-simple-problem)
   - [6.2 Adaptive Dormand-Prince on a hard problem (Lorenz)](#62-adaptive-dormand-prince-on-a-hard-problem)
   - [6.3 Event detection — artillery shot terminated on impact](#63-event-detection--artillery-shot)
   - [6.4 Boundary value problem — shooting method](#64-boundary-value-problem--shooting-method)
   - [6.5 Where to aim? — shooting with a free boundary](#65-where-to-aim--shooting-with-a-free-boundary)
   - [6.6 Stiff system — Rosenbrock 2(3)](#66-stiff-system--rosenbrock-23)
7. ✅ **[Eigenvalues and eigenvectors](#recipe-7-eigenvalues-and-eigenvectors)**
   - [7.1 Eigenvalues and eigenvectors of a real symmetric matrix](#71-real-symmetric-matrix)
   - [7.2 Eigenvalues of a general real matrix](#72-general-real-matrix)
   - [7.3 Eigenvalues and eigenvectors of a complex Hermitian matrix](#73-complex-hermitian-matrix)
   - [7.4 Eigenvalues of a general complex matrix](#74-general-complex-matrix)
8. ✅ **[Matrix decompositions](#recipe-8-matrix-decompositions)**
    - [8.1 SVD: singular values, rank, and condition number](#81-svd-singular-values-rank-and-condition-number)
    - [8.2 SVD: low-rank approximation](#82-svd-low-rank-approximation)
    - [8.3 SVD: null space and rank deficiency](#83-svd-null-space-and-rank-deficiency)
    - [8.4 QR: orthogonal basis and reconstruction](#84-qr-orthogonal-basis-and-reconstruction)
    - [8.5 LU: factors and determinant](#85-lu-factors-and-determinant)
    - [8.6 Cholesky: positive-definite decomposition](#86-cholesky-positive-definite-decomposition)
9. ✅ **[Matrix properties and diagnostics](#recipe-9-matrix-properties-and-diagnostics)**
    - [9.1 Inspect matrix structure](#91-inspect-matrix-structure)
    - [9.2 Norms and approximation error](#92-norms-and-approximation-error)
    - [9.3 Trace, determinant, and rank](#93-trace-determinant-and-rank)
    - [9.4 Condition number and sensitivity](#94-condition-number-and-sensitivity)
    - [9.5 Inverse and identity verification](#95-inverse-and-identity-verification)
    - [9.6 Matrix slicing and transformations](#96-matrix-slicing-and-transformations)
    - [9.7 Matrix-valued functions](#97-matrix-valued-functions)
10. ✅ **[Statistics](#recipe-10-statistics)**
    - [10.1 Quick descriptive summary](#101-quick-descriptive-summary)
    - [10.2 Robust statistics with outliers](#102-robust-statistics-with-outliers)
    - [10.3 Weighted measurements](#103-weighted-measurements)
    - [10.4 Relationships between variables](#104-relationships-between-variables)
    - [10.5 Histogram and empirical CDF](#105-histogram-and-empirical-cdf)
    - [10.6 Confidence interval for an unknown mean](#106-confidence-interval-for-an-unknown-mean)
    - [10.7 Real dataset — Palmer Penguins](#107-real-dataset--palmer-penguins)
    - [10.8 Real analysis — tyre degradation at Monza](#108-real-analysis--tyre-degradation-at-monza)
11. ✅ **[3D geometry: planes, lines, and bodies](#recipe-11-3d-geometry-planes-lines-and-bodies)**
    - [11.1 Build points, vectors, lines, and segments](#111-build-points-vectors-lines-and-segments)
    - [11.2 Intersect a line with a plane](#112-intersect-a-line-with-a-plane)
    - [11.3 Project points and measure distances](#113-project-points-and-measure-distances)
    - [11.4 Work with triangles and rectangular surfaces](#114-work-with-triangles-and-rectangular-surfaces)
    - [11.5 Create and query 3D bodies](#115-create-and-query-3d-bodies)
    - [11.6 Use bounding volumes for spatial checks](#116-use-bounding-volumes-for-spatial-checks)
12. ✅ **[3D geometry: algorithms](#recipe-12-3d-geometry-algorithms)**
    - [12.1 Classify two 3D lines](#121-classify-two-3d-lines)
    - [12.2 Find closest approach between skew lines](#122-find-closest-approach-between-skew-lines)
    - [12.3 Project points onto lines, segments, and planes](#123-project-points-onto-lines-segments-and-planes)
    - [12.4 Intersect lines with planes and planes with planes](#124-intersect-lines-with-planes-and-planes-with-planes)
    - [12.5 Work with triangle geometry and point containment](#125-work-with-triangle-geometry-and-point-containment)
    - [12.6 Ray-pick a triangle](#126-ray-pick-a-triangle)
    - [12.7 Query a 3D point cloud efficiently](#127-query-a-3d-point-cloud-efficiently)
    - [12.8 Build a convex hull and test containment](#128-build-a-convex-hull-and-test-containment)
13. ✅ **[Interpolating functions](#recipe-13-interpolating-functions)**
    - [13.1 Interpolate a real function with several methods](#131-interpolate-a-real-function-with-several-methods)
    - [13.2 Interpolate a 2D function on a grid](#132-interpolate-a-2d-function-on-a-grid)
    - [13.3 Interpolate a parametric curve](#133-interpolate-a-parametric-curve)
14. ✅ **[Curve fitting, including fitting to points](#recipe-14-curve-fitting-including-fitting-to-points)**
    - [14.1 Fit a sensor calibration line](#141-fit-a-sensor-calibration-line)
    - [14.2 Use weighted least squares for unequal measurement uncertainty](#142-use-weighted-least-squares-for-unequal-measurement-uncertainty)
    - [14.3 Fit a polynomial performance curve](#143-fit-a-polynomial-performance-curve)
    - [14.4 Fit seasonal demand with Fourier basis functions](#144-fit-seasonal-demand-with-fourier-basis-functions)
    - [14.5 Fit an exponential response curve](#145-fit-an-exponential-response-curve)
    - [14.6 Fit a smooth 2D path from points](#146-fit-a-smooth-2d-path-from-points)
15. ✅ **[Curves and surfaces](#recipe-15-curves-and-surfaces)**
    - [15.1 Use built-in curves and surfaces](#151-use-built-in-curves-and-surfaces)
    - [15.2 Define custom curves and surfaces](#152-define-custom-curves-and-surfaces)
    - [15.3 Calculate basic curve properties](#153-calculate-basic-curve-properties)
    - [15.4 Compute a Frenet frame](#154-compute-a-frenet-frame)
    - [15.5 Inspect surface normals and curvatures](#155-inspect-surface-normals-and-curvatures)
    - [15.6 Calculate first and second fundamental forms](#156-calculate-first-and-second-fundamental-forms)
  16. ✅ **[Polynoms: polynomial models and approximations](#recipe-16-polynoms-polynomial-models-and-approximations)**
      - [16.1 Build and evaluate a polynomial model](#161-build-and-evaluate-a-polynomial-model)
      - [16.2 Compose polynomial arithmetic into a power model](#162-compose-polynomial-arithmetic-into-a-power-model)
      - [16.3 Divide polynomials to inspect quotient and residual](#163-divide-polynomials-to-inspect-quotient-and-residual)
      - [16.4 Differentiate and integrate polynomial models](#164-differentiate-and-integrate-polynomial-models)
      - [16.5 Construct an explicit polynomial from measurements](#165-construct-an-explicit-polynomial-from-measurements)
      - [16.6 Use Chebyshev approximation and convert to a polynomial](#166-use-chebyshev-approximation-and-convert-to-a-polynomial)
17. ✅ **[Quaternions](#recipe-17-quaternions)**
    - [17.1 Rotate a vector around an arbitrary axis](#171-rotate-a-vector-around-an-arbitrary-axis)
    - [17.2 Point one direction toward another](#172-point-one-direction-toward-another)
    - [17.3 Compose rotations in the correct order](#173-compose-rotations-in-the-correct-order)
    - [17.4 Convert between Euler angles, quaternions, and matrices](#174-convert-between-euler-angles-quaternions-and-matrices)
    - [17.5 Undo rotations and calculate relative orientation](#175-undo-rotations-and-calculate-relative-orientation)
    - [17.6 Interpolate orientations with SLERP](#176-interpolate-orientations-with-slerp)
18. ✅ **[Analyzing functions](#recipe-18-analyzing-functions)**
    - [18.1 Produce a point and interval health report](#181-produce-a-point-and-interval-health-report)
    - [18.2 Build a complete polynomial feature map](#182-build-a-complete-polynomial-feature-map)
    - [18.3 Detect and classify discontinuities](#183-detect-and-classify-discontinuities)
    - [18.4 Estimate oscillation period from zero crossings](#184-estimate-oscillation-period-from-zero-crossings)
    - [18.5 Measure approximation error](#185-measure-approximation-error)
19. ✅ **[Coordinate transformations](#recipe-19-coordinate-transformations)** · *“coordinates → curvature” thread, part 1 of 4*
    - [19.1 Points: Cartesian ↔ spherical ↔ cylindrical](#191-points-cartesian--spherical--cylindrical)
    - [19.2 Vectors are not points — contravariant velocity](#192-vectors-are-not-points-a-velocity-transforms-contravariantly)
    - [19.3 Gradients transform covariantly](#193-gradients-transform-by-the-other-rule-covariantly)
20. ✅ **[Fields](#recipe-20-fields)** · *thread part 2 of 4*
    - [20.1 One field, two coordinate systems, one answer](#201-one-field-two-coordinate-systems-one-physical-answer)
    - [20.2 curl grad = 0, div curl = 0](#202-the-identities-every-field-obeys)
    - [20.3 ∇²(1/r) = 0 — harmonic potential](#203-empty-space-gravity-1r--0)
21. ✅ **[Tensors](#recipe-21-tensors)** · *thread part 3 of 4*
    - [21.1 Tensor2 basics: variance, contraction, evaluation](#211-tensor2-basics-index-variance-contraction-evaluation)
    - [21.2 The metric — what makes components meaningful](#212-the-tensor-the-metric--what-makes-components-meaningful)
    - [21.3 Raising an index: gradient → force](#213-raising-an-index-gradient-covariant--force-contravariant)
    - [21.4 The twist: curvy coordinates, flat space](#214-the-twist-curvy-coordinates-flat-space)
22. ✅ **[Differential geometry and geodesics](#recipe-22-differential-geometry-and-geodesics)** · *thread finale*
    - [22.1 Induced metric + Theorema Egregium](#221-the-spheres-induced-metric-and-the-theorema-egregium)
    - [22.2 Geodesics: the great circles](#222-geodesics-the-great-circles)
    - [22.3 Finale — holonomy: curvature made visible](#223-finale--holonomy-parallel-transport-remembers-curvature)

### Advanced Recipes

23. ✅ **[Fourier algorithms, including FFT](#recipe-23-fourier-algorithms)**
    - [23.1 FFT round trip and normalization](#231-fft-round-trip-and-normalization)
    - [23.2 Detect frequencies and amplitudes](#232-detect-frequencies-and-amplitudes)
    - [23.3 Control spectral leakage with windows](#233-control-spectral-leakage-with-windows)
    - [23.4 Remove interference in the frequency domain](#234-remove-interference-in-the-frequency-domain)
    - [23.5 Smooth signals with convolution](#235-smooth-signals-with-convolution)
    - [23.6 Estimate time delay from phase](#236-estimate-time-delay-from-phase)
    - [23.7 Transform arbitrary sample counts with Bluestein](#237-transform-arbitrary-sample-counts-with-bluestein)
    - [23.8 Verify energy with Parseval's theorem](#238-verify-energy-with-parsevals-theorem)
    - [23.9 Compress smooth data with the DCT](#239-compress-smooth-data-with-the-dct)
    - [23.10 Analyze zero-boundary modes with the DST](#2310-analyze-zero-boundary-modes-with-the-dst)
24. ✅ **[Graph algorithms](#recipe-24-graph-algorithms)**
    - [24.1 Build, traverse, and inspect a graph](#241-build-traverse-and-inspect-a-graph)
    - [24.2 Find the best route](#242-find-the-best-route)
    - [24.3 Schedule dependent tasks and find the critical path](#243-schedule-dependent-tasks-and-find-the-critical-path)
    - [24.4 Design the cheapest connected network](#244-design-the-cheapest-connected-network)
    - [24.5 Find single points of failure](#245-find-single-points-of-failure)
    - [24.6 Analyze cycles and strongly connected subsystems](#246-analyze-cycles-and-strongly-connected-subsystems)
    - [24.7 Calculate maximum throughput and the bottleneck cut](#247-calculate-maximum-throughput-and-the-bottleneck-cut)
    - [24.8 Assign jobs to workers](#248-assign-jobs-to-workers)
25. ✅ **[DAE solvers](#recipe-25-dae-solvers)**
    - [25.1 Define and solve an index-1 DAE](#251-define-and-solve-an-index-1-dae)
    - [25.2 Compute consistent initial algebraic variables](#252-compute-consistent-initial-algebraic-variables)
    - [25.3 Compare Backward Euler, BDF2/BDF4, and Radau IIA](#253-compare-backward-euler-bdf2bdf4-and-radau-iia)
    - [25.4 Handle DAE events and restart with consistent constraints](#254-handle-dae-events-and-restart-with-consistent-constraints)
26. ✅ **[Algebra: finite structures, symmetry, and exact arithmetic](#recipe-26-algebra-finite-structures-symmetry-and-exact-arithmetic)**
    - [26.1 Reorder real data with permutations](#261-reorder-real-data-with-permutations)
    - [26.2 Model polygon symmetries with a dihedral group](#262-model-polygon-symmetries-with-a-dihedral-group)
    - [26.3 Count distinct bracelets with Burnside's lemma](#263-count-distinct-bracelets-with-burnsides-lemma)
    - [26.4 Work exactly with modular arithmetic and prime fields](#264-work-exactly-with-modular-arithmetic-and-prime-fields)
    - [26.5 Solve an exact linear system over a finite field](#265-solve-an-exact-linear-system-over-a-finite-field)
    - [26.6 Build a tiny extension field from polynomials](#266-build-a-tiny-extension-field-from-polynomials)
27. ✅ **[Complex analysis](#recipe-27-complex-analysis)**
    - [27.1 Differentiate an analytic complex function](#271-differentiate-an-analytic-complex-function)
    - [27.2 Integrate along complex contours](#272-integrate-along-complex-contours)
    - [27.3 Compute winding numbers](#273-compute-winding-numbers)
    - [27.4 Recover values and derivatives with Cauchy's formula](#274-recover-values-and-derivatives-with-cauchys-formula)
    - [27.5 Compute residues and verify the residue theorem](#275-compute-residues-and-verify-the-residue-theorem)
    - [27.6 Count zeros and poles with the argument principle](#276-count-zeros-and-poles-with-the-argument-principle)
    - [27.7 Count roots first, then locate them](#277-count-roots-first-then-locate-them)
28. Function spaces *(planned)*
29. Vector spaces *(planned)*
30. ✅ **[Combinatorics and number theory](#recipe-30-combinatorics-and-number-theory)**
    - [30.1 Select a committee and assign roles](#301-select-a-committee-and-assign-roles)
    - [30.2 Enumerate constrained resource allocations](#302-enumerate-constrained-resource-allocations)
    - [30.3 Count ways to group labeled objects](#303-count-ways-to-group-labeled-objects)
    - [30.4 Analyze an integer through its prime structure](#304-analyze-an-integer-through-its-prime-structure)
    - [30.5 Educational RSA-style encryption](#305-educational-rsa-style-encryption)
    - [30.6 Synchronize repeating schedules with CRT](#306-synchronize-repeating-schedules-with-crt)
31. Using special functions *(planned)*
32. Serialization *(planned)*
33.  Persistence *(planned)*
34. Visualization *(planned)*i

---

## Recipe 1: Solving linear systems

### 1.1 Real 3×3 system

**Problem:** solve the 3×3 system Ax = b

```
 2x +  y −  z =   8
−3x −  y + 2z = −11
−2x +  y + 2z =  −3
```

(the classic textbook system; exact solution x = (2, 3, −1)).

**Setup**

```cpp
#include <mml/core/LinAlgEqSolvers/LinAlgDirect.h>
#include <mml/systems/LinearSystem.h>

using namespace MML;

Matrix<Real> A(3, 3, {  2,  1, -1,
                       -3, -1,  2,
                       -2,  1,  2 });
Vector<Real> b({ 8, -11, -3 });
```

**Solve**

Directly, with an explicit solver — `LUSolver` decomposes once (partial pivoting) and can
then solve for any number of right-hand sides:

```cpp
LUSolver<Real> lu(A);              // decomposition happens here; throws if singular
Vector<Real> x = lu.Solve(b);      // x = (2, 3, -1)
```

With automatic solver selection — the one-liner:

```cpp
Vector<Real> x2 = Systems::SolveLinearSystem(A, b);
```

Or the `LinearSystem` facade when you want diagnostics and control:

```cpp
Systems::LinearSystem<Real> sys(A, b);

Vector<Real> x3 = sys.Solve();                 // auto-selects the best solver
auto check = sys.Verify(x3);                   // residual-based quality check
// check.isAccurate == true

Vector<Real> x4 = sys.SolveByQR();             // or force a specific method
```

### 1.2 Complex 3×3 system

The direct solvers are templates over the element type, so a complex system is the *same*
recipe with `Complex` elements — no separate API to learn.

**Setup** — we manufacture b from a known solution x = (1, i, 1−i), so the answer is
verifiable by construction:

```cpp
Matrix<Complex> Ac(3, 3, { Complex(2, 1), Complex(1, 0), Complex(0, -1),
                           Complex(0, 2), Complex(3, -1), Complex(1, 0),
                           Complex(1, 0), Complex(0, 1), Complex(2, 2) });

Vector<Complex> xExpected({ Complex(1, 0), Complex(0, 1), Complex(1, -1) });
Vector<Complex> bc = Ac * xExpected;
```

**Solve**

```cpp
LUSolver<Complex> luc(Ac);
Vector<Complex> xc = luc.Solve(bc);
// xc.IsEqualTo(xExpected, 1e-12) == true
```

**Notes**

- `LUSolver` = decompose once, reuse for many right-hand sides; `GaussJordanSolver` and
  `CholeskySolver` (SPD matrices) follow the same shape.
- `Systems::LinearSystem::Solve()` analyzes the matrix and picks the appropriate
  decomposition (triangular / Cholesky / QR / SVD / LU); `SolveBy...()` methods force one.
  The facade is `Real`-only — complex systems use the direct solvers as in 1.2.
- `Verify(x)` returns a `VerificationResult` with absolute/relative residuals, backward
  error, and the `isAccurate` verdict.
- Singular or inconsistent systems throw MML exceptions from the direct solvers;
  use `sys.Analyze()` first when the system might be degenerate.
- Runnable version: `Cookbook_Recipe01_SolvingLinearSystems()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 2: Root finding

All functions here come from `MML::RootFinding` (aggregate header
`<mml/algorithms/RootFinding.h>`).

### 2.1 Classic: bracket the root, then solve

**Problem:** solve **Kepler's equation** E − e·sin E = M for the eccentric anomaly E of a
highly eccentric orbit (e = 0.9, M = 1). You rarely know a bracketing interval up front —
so let MML find one first.

**Setup**

```cpp
#include <mml/algorithms/RootFinding.h>

RealFunction kepler([](Real E) { return E - 0.9 * std::sin(E) - 1.0; });
```

**Solve**

```cpp
// Expand an initial interval outward until it brackets a sign change...
Real x1 = 0.0, x2 = 0.5;
bool bracketed = RootFinding::BracketRoot(kepler, x1, x2);   // x1, x2 now bracket the root

// ...then hand the bracket to a robust solver
RootFinding::RootFindingConfig config;
config.tolerance = 1e-14;
auto result = RootFinding::FindRootBrent(kepler, x1, x2, config);
// result.root = 1.86208... (E for e=0.9, M=1), result.converged, result.iterations_used
```

`FindRootBrackets` is the sibling that scans an interval and returns *all* sub-brackets —
the manual road to §2.2.

### 2.2 All roots in an interval

**Problem:** a vibrating cantilever beam has its natural frequencies where
cos(x)·cosh(x) + 1 = 0 — an oscillatory equation with many roots.

```cpp
RealFunction beam([](Real x) { return std::cos(x) * std::cosh(x) + 1.0; });

auto all = RootFinding::FindAllRealRootsInInterval(beam, 0.0, 12.0);
// all.roots = { 1.87510, 4.69409, 7.85476, 10.99554 }  - the classic beam eigenvalues
```

The result also carries per-root refinement diagnostics (`all.accepted`) and rejected
candidates; tune isolation density via `FindAllRealRootsConfig::isolation.num_intervals`
(default 100 samples).

### 2.3 Polynomial roots

**Problem:** all roots of p(x) = (x−1)(x−2)(x−3)(x−4)(x−5) — the “Wilkinson-5”,
notoriously sensitive to coefficient errors.

```cpp
#include <mml/base/Polynom.h>
#include <mml/algorithms/RootFinding/RootFindingPolynoms.h>

// Coefficients in ASCENDING power order: a0 + a1*x + ... + a5*x^5
PolynomReal p({ -120, 274, -225, 85, -15, 1 });

Vector<Complex> roots = RootFinding::LaguerreRoots(p);   // 1, 2, 3, 4, 5
```

Polynomial roots are complex in general, so the result is `Vector<Complex>` even for a
real polynomial. Alternatives: `EigenvalueRoots` (companion-matrix) and `BairstowRoots`.

### 2.4 Roots of a complex function

**Problem:** eᶻ = z has **no real solutions** — its roots live properly in ℂ
(the principal one is z ≈ 0.31813 + 1.33724i).

```cpp
#include <mml/base/ComplexFunction.h>
#include <mml/algorithms/RootFinding/RootFindingComplex.h>

ComplexFunction f([](Complex z) { return std::exp(z) - z; });

// Muller's method from a single complex starting point
auto result = RootFinding::FindRootMuller(f, Complex(0.5, 1.0));
// result.root = (0.318131, 1.337236), result.converged, |result.function_value| ~ 1e-11
```

`FindRootNewtonComplex` is the Newton alternative (numerical complex derivative).

### 2.5 Roots in multiple dimensions

**Problem:** find an intersection of the circle x² + y² = 4 with the curve eˣ + y = 1 —
a 2×2 nonlinear system F(x, y) = 0.

```cpp
RootFinding::DynamicSystemFunction F = [](const Vector<Real>& p) {
    return Vector<Real>({ p[0] * p[0] + p[1] * p[1] - 4.0,
                          std::exp(p[0]) + p[1] - 1.0 });
};

// Damped Newton; the Jacobian is computed numerically when not supplied
auto result = RootFinding::SolveNonlinearSystemNewton(F, Vector<Real>({ 1.0, -1.7 }));
// result.solution = (1.004169, -1.729637), result.residual_norm, result.IsSuccess()
```

Supply an analytic Jacobian (`DynamicJacobianFunction`) for speed; fixed-size variants
(`VectorN<Real, N>` + `IVectorFunction<N>` overloads) avoid heap allocation, and
`config.store_trace = true` records the full iteration history.

**Notes**

- Every solver returns the standard Result object (Fundamentals §4): check `converged`
  (scalar/complex) or `IsSuccess()` (systems) before using the answer.
- Runnable version: `Cookbook_Recipe02_RootFinding()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 3: Calculating integrals

Everything below comes from the aggregate header `<mml/core/Integration.h>` unless noted;
path/surface/Monte-Carlo live in their own headers under `mml/core/Integration/`.

### 3.1 Real function — trapezoid, Romberg, Gauss-Kronrod

**Problem:** ∫₀^{2π} e^{−x/2} cos 3x dx — a damped oscillation with a known antiderivative,
so every method can be checked against the exact value (1 − e^{−π})/2 / 9.25 ≈ 0.0517182.

```cpp
RealFunction f([](Real x) { return std::exp(-x / 2) * std::cos(3 * x); });
const Real a = 0.0, b = 2 * Constants::PI;

Real trap    = IntegrateTrap(f, a, b);       // workhorse; IntegrationResult -> Real
Real romberg = IntegrateRomberg(f, a, b);    // extrapolated - fast for smooth f

Real gkError = 0.0;                          // adaptive Gauss-Kronrod + error estimate
Real gk = Integration::IntegrateGaussKronrod(f, a, b, &gkError);
```

Measured accuracy on this integrand: trapezoid ~1e-6, Romberg ~5e-11, Gauss-Kronrod exact
to machine precision with a trustworthy error estimate (`gkError` ≈ 4e-11). There are also
`IntegrateSimpson`, `IntegrateGauss10`, `IntegrateGK21`, and a default `Integrate(f, a, b, eps)`.

### 3.2 2D integration — simple and adaptive

**Simple** handles curved domains — the y-limits are functions of x. Integral of x²+y²
over the unit disk (= π/2):

```cpp
ScalarFunction<2> f2([](const VectorN<Real, 2>& p) { return p[0] * p[0] + p[1] * p[1]; });

Real disk = Integrate2D(f2, IntegrationMethod::GAUSS10,
    -1, 1,
    [](Real x) { return -std::sqrt(1 - x * x); },
    [](Real x) { return  std::sqrt(1 - x * x); });
```

**Adaptive** works on rectangles with error control and diagnostics — it takes a plain
`(Real x, Real y)` callable; ∫∫ sin x · cos y over [0,π]×[0,π/2] (= 2):

```cpp
auto g2 = [](Real x, Real y) { return std::sin(x) * std::cos(y); };

auto adaptive = Integration::IntegrateAdaptive2D(g2, 0.0, Constants::PI, 0.0, Constants::PI / 2, 1e-10);
// adaptive.value == 2 to 1e-12, adaptive.error_estimate, adaptive.function_evaluations (225)
```

### 3.3 3D integration — simple and adaptive

Same pattern one dimension up. Volume of the unit sphere via nested limit lambdas:

```cpp
ScalarFunction<3> one([](const VectorN<Real, 3>&) { return Real(1); });

Real vol = Integrate3D(one,
    -1, 1,
    [](Real x) { return -std::sqrt(1 - x * x); },
    [](Real x) { return  std::sqrt(1 - x * x); },
    [](Real x, Real y) { return -std::sqrt(std::max(Real(0), 1 - x * x - y * y)); },
    [](Real x, Real y) { return  std::sqrt(std::max(Real(0), 1 - x * x - y * y)); });
// vol = 4.1907  (exact 4pi/3 = 4.18879)
```

Adaptive over a box, `(Real x, Real y, Real z)` callable; ∫∫∫ xyz over [0,1]³ (= 1/8):

```cpp
auto g3 = [](Real x, Real y, Real z) { return x * y * z; };
auto adaptive = Integration::IntegrateAdaptive3D(g3, 0, 1, 0, 1, 0, 1, 1e-10);   // 0.125 exactly
```

### 3.4 Improper integrals

Three classics — infinite upper limit, whole real line, and an integrable endpoint
singularity:

```cpp
// Gaussian tail:  ∫ 0..inf e^(-x^2) dx = sqrt(pi)/2
RealFunction gaussian([](Real x) { return std::exp(-x * x); });
Real tail = IntegrateUpperInf(gaussian, 0.0);                    // 0.886227

// Lorentzian:  ∫ -inf..inf 1/(1+x^2) dx = pi
RealFunction lorentzian([](Real x) { return 1.0 / (1.0 + x * x); });
Real whole = IntegrateInf(lorentzian);                           // 3.141593

// Endpoint singularity:  ∫ 0..1 1/sqrt(x) dx = 2
RealFunction invSqrt([](Real x) { return 1.0 / std::sqrt(x); });
Real singular = IntegrateLowerSingular(invSqrt, 0.0, 1.0);       // 2
```

Siblings: `IntegrateLowerInf`, `IntegrateInfSplit`, `IntegrateUpperSingular`,
`IntegrateBothSingular`, and open-interval `IntegrateOpen`.

### 3.5 Path integration — work done in a force field

**Problem:** work W = ∫ F·dr done by the force F = (−y, x, z) along one turn of the helix
(cos t, sin t, t). Exact: W = 2π + 2π² ≈ 26.0224.

```cpp
#include <mml/core/Integration/PathIntegration.h>

VectorFunction<3> force([](const VectorN<Real, 3>& p) {
    return VectorN<Real, 3>{ -p[1], p[0], p[2] };
});
ParametricCurve<3> helix([](Real t) { return VectorN<Real, 3>{ std::cos(t), std::sin(t), t }; });

Real work = PathIntegration::LineIntegral(force, helix, 0.0, 2 * Constants::PI, 1e-8);
// work = 26.02239411 - matches 2pi + 2pi^2 to all displayed digits
```

Scalar line integrals (∫ f ds), `ParametricCurveLength`, and `ParametricCurveMass` live in
the same class.

### 3.6 Surface integration — flux through a closed surface

**Problem:** flux Φ = ∮ F·dS of the radial field F = (x, y, z) through a cube of side 2
centered at the origin. By the divergence theorem Φ = ∫ div F dV = 3·V = 24.

```cpp
#include <mml/base/Geometry/Geometry3DBodies.h>
#include <mml/core/Integration/SurfaceIntegration.h>

VectorFunction<3> field([](const VectorN<Real, 3>& p) { return p; });
Cube3D cube(2.0);

Real flux = SurfaceIntegration::SurfaceIntegral(field, cube);    // 24, exactly
```

The same call accepts any `BodyWithRectSurfaces` / `BodyWithTriangleSurfaces` solid or a
single parametric surface — this is the machinery behind the README's Gauss-divergence
flagship example.

### 3.7 Monte Carlo integration

**Problem:** ∫∫∫ e^{−|p|²} over [−1,1]³; exact value (√π · erf 1)³ ≈ 3.33231.

```cpp
#include <mml/core/Integration/MonteCarloIntegration.h>

ScalarFunction<3> gauss3([](const VectorN<Real, 3>& p) {
    return std::exp(-(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]));
});

MonteCarloIntegrator<3> mc(42);                 // seeded - reproducible
MonteCarloConfig config;
config.num_samples = 500000;

auto result = mc.integrate(gauss3,
    VectorN<Real, 3>{ -1, -1, -1 }, VectorN<Real, 3>{ 1, 1, 1 }, config);
// result.value = 3.33203 (within one std error 0.0023 of exact), result.samples_used
```

`StratifiedMonteCarloIntegrator` and `HitOrMissIntegrator` are the variance-reduction
siblings; `EstimatePi()` and `EstimateUnitBallVolume()` are ready-made classics.

**Notes**

- 1D quadratures return `IntegrationResult` (implicitly convertible to `Real`); adaptive
  and Monte Carlo results carry error estimates and evaluation counts — read them.
- Improper integrals rely on variable transforms internally; pick the routine that matches
  where the trouble is (infinity vs. endpoint singularity).
- Runnable version: `Cookbook_Recipe03_Integrals()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 4: Calculating derivatives

The numerical differentiation API lives in `<mml/core/Derivation.h>`. One naming detail is
essential: in `NDer1`, `NDer2`, `NDer4`, `NDer6`, and `NDer8`, the number describes the
finite-difference **accuracy order**, not the derivative order. Actual second and third
derivatives use the `NSecDer...` and `NThirdDer...` families.

### 4.1 Real function — first, second, and third derivatives

**Problem:** calculate the first three derivatives of $f(x)=\sin x$ at $x=\pi/4$. The exact
values are $f'=\sqrt{2}/2$, $f''=-\sqrt{2}/2$, and $f'''=-\sqrt{2}/2$.

**Setup**

```cpp
#include <mml/core/Derivation.h>

RealFunction f([](Real x) { return std::sin(x); });
const Real x = Constants::PI / 4.0;
```

**Solve** — use fourth-order-accurate formulas for each derivative order:

```cpp
Real first  = Derivation::NDer4(f, x);
Real second = Derivation::NSecDer4(f, x);
Real third  = Derivation::NThirdDer4(f, x);
// first = 0.70710678, second = -0.70710678, third = -0.70710684
```

The overloads choose scale-aware step sizes automatically. Pass an explicit `h` before the
optional error pointer when the problem has a known numerical scale.

### 4.2 Scalar function — one partial and all partials

**Problem:** for $f(x,y,z)=x^2y+\sin z$, calculate $\partial f/\partial x$ and the full
gradient at $(1,2,0)$. Analytically, $\nabla f=(2xy,x^2,\cos z)=(4,1,1)$.

**Setup**

```cpp
ScalarFunction<3> f([](const VectorN<Real, 3>& p) {
    return p[0] * p[0] * p[1] + std::sin(p[2]);
});
VectorN<Real, 3> point{ 1.0, 2.0, 0.0 };
```

**Solve**

```cpp
Real df_dx = Derivation::NDer4Partial(f, 0, point);
VectorN<Real, 3> gradient = Derivation::NDer4PartialByAll(f, point);
// df_dx = 4, gradient = (4, 1, 1)
```

Indices are zero-based. `NDer4Partial` differentiates with respect to one selected input;
`NDer4PartialByAll` returns all first partials, which is the gradient of a scalar function.

### 4.3 Vector function — partial, partial-by-all, all-by-all

**Problem:** differentiate $F(x,y,z)=(xy,yz,zx)$ at $(1,2,3)$. A vector function has two
indices: which output component $F_i$, and which input variable $x_j$.

**Setup**

```cpp
VectorFunction<3> field([](const VectorN<Real, 3>& p) {
    return VectorN<Real, 3>{ p[0] * p[1], p[1] * p[2], p[2] * p[0] };
});
VectorN<Real, 3> point{ 1.0, 2.0, 3.0 };
```

**Solve**

```cpp
Real partial = Derivation::NDer4Partial(field, 0, 1, point); // dF0/dy = 1

VectorN<Real, 3> firstComponentGradient =
    Derivation::NDer4PartialByAll(field, 0, point);           // (2, 1, 0)

MatrixNM<Real, 3, 3> allByAll =
    Derivation::NDer4PartialAllByAll(field, point);           // full Jacobian
```

`Partial` returns one $\partial F_i/\partial x_j$, `PartialByAll` differentiates one output
with respect to every input, and `PartialAllByAll` returns every pair as a matrix whose rows
are output components and columns are input variables.

### 4.4 Parametric curve — first, second, and third derivatives

**Problem:** calculate velocity, acceleration, and jerk along the circular helix
$r(t)=(\cos t,\sin t,t)$ at $t=0$.

**Setup**

```cpp
ParametricCurve<3> helix([](Real t) {
    return VectorN<Real, 3>{ std::cos(t), std::sin(t), t };
});
```

**Solve**

```cpp
VectorN<Real, 3> first  = Derivation::NDer4(helix, 0.0);
VectorN<Real, 3> second = Derivation::NSecDer4(helix, 0.0);
VectorN<Real, 3> third  = Derivation::NThirdDer4(helix, 0.0);
// r'(0) = (0, 1, 1), r''(0) = (-1, 0, 0), r'''(0) = (0, -1, 0)
```

These return vectors directly. The second- and third-derivative routines use direct finite
difference formulas rather than repeatedly differentiating an already approximate result.

### 4.5 Parametric surface — tangent plane, normal, and area scale

**Problem:** at one point on a sphere of radius 2, numerically recover the tangent plane,
outward unit normal, and the local area scale $\lVert r_u\times r_w\rVert$.

**Setup** — $u$ is polar angle and $w$ is azimuth:

```cpp
ParametricSurfaceRect<3> sphere([](Real u, Real w) {
    return VectorN<Real, 3>{
        2.0 * std::sin(u) * std::cos(w),
        2.0 * std::sin(u) * std::sin(w),
        2.0 * std::cos(u)
    };
});
const Real u = Constants::PI / 4.0;
const Real w = Constants::PI / 3.0;
```

**Solve** — the two parameter derivatives span the tangent plane; their cross product is
normal to it and its magnitude converts parameter-space area $du\,dw$ to surface area:

```cpp
VectorN<Real, 3> tangentU = Derivation::NDer2_u(sphere, u, w);
VectorN<Real, 3> tangentW = Derivation::NDer2_w(sphere, u, w);
VectorN<Real, 3> normalCross{
    tangentU[1] * tangentW[2] - tangentU[2] * tangentW[1],
    tangentU[2] * tangentW[0] - tangentU[0] * tangentW[2],
    tangentU[0] * tangentW[1] - tangentU[1] * tangentW[0]
};
Real areaScale = normalCross.NormL2();
VectorN<Real, 3> unitNormal = normalCross / areaScale;
// unitNormal = (0.353553, 0.612372, 0.707107), areaScale = 2.828427
```

Here `NDer2_u/w` means a second-order-accurate **first** derivative. For curvature work,
`NDer2_uu`, `NDer2_uw`, and `NDer2_ww` provide the true second parameter derivatives.

### 4.6 Jacobians and Hessians

**Problem:** build the complete Jacobian of $F(x,y)=(x^2+y,xy)$ and the Hessian of
$f(x,y)=x^2+xy+3y^2$ at $(1,2)$.

**Setup**

```cpp
VectorFunction<2> map([](const VectorN<Real, 2>& p) {
    return VectorN<Real, 2>{ p[0] * p[0] + p[1], p[0] * p[1] };
});
ScalarFunction<2> quadratic([](const VectorN<Real, 2>& p) {
    return p[0] * p[0] + p[0] * p[1] + 3.0 * p[1] * p[1];
});
VectorN<Real, 2> point{ 1.0, 2.0 };
```

**Solve**

```cpp
MatrixNM<Real, 2, 2> jacobian = Derivation::calcJacobian(map, point);
MatrixNM<Real, 2, 2> hessian = Derivation::calcHessian(quadratic, point);
// Jacobian = [[2, 1], [2, 1]]; Hessian = [[2, 1], [1, 6]]
```

The Jacobian convention is $J_{ij}=\partial F_i/\partial x_j$. The Hessian is assembled
symmetrically from fourth-order formulas for both pure and mixed second partials.

### 4.7 High-precision complex-step differentiation

**Problem:** differentiate $f(x)=e^x\sin x$ at $x=1$ using an extremely small step. An
ordinary real finite difference subtracts nearly equal numbers and collapses to zero;
complex-step differentiation avoids that subtraction.

**Setup** — write one analytic formula that accepts either `Real` or `Complex`:

```cpp
auto f = [](auto x) { return std::exp(x) * std::sin(x); };
const Real x = 1.0;
const Real exact = std::exp(x) * (std::sin(x) + std::cos(x));
```

**Solve**

```cpp
Real tinyFiniteDifference = Derivation::NDer2(f, x, 1e-20);
Real complexStep = Derivation::ComplexStep(f, x); // default h is extremely small
// exact = 3.756049..., tinyFiniteDifference = 0, complexStep = 3.756049...
```

The identity $f'(x)=\operatorname{Im}(f(x+ih))/h+O(h^2)$ has no subtractive cancellation,
so `h` can be far below the useful finite-difference range. It requires an analytic
continuation implemented with complex-aware operations; non-analytic operations such as
`abs`, `real`, branching on the argument, or conjugation generally invalidate the method.

**Notes**

- The default derivative step is scale-aware; an arbitrarily smaller finite-difference step
  is usually less accurate because roundoff eventually dominates truncation error.
- Pass an error pointer to the finite-difference routines when an estimate is useful, or use
  the corresponding `...Detailed` real-function API for status and diagnostics.
- `NDer4PartialAllByAll` and `calcJacobian` compute the same mathematical object; the latter
  is the clearer intent-level API when the whole Jacobian is wanted.
- Runnable version: `Cookbook_Recipe04_Derivatives()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 5: Optimization recipes

MML separates scalar one-dimensional minimization, multidimensional unconstrained
minimization, and dense linear programming. Use `<mml/algorithms/Optimization/Optimization.h>`
for golden section and Brent, `<mml/algorithms/Optimization/OptimizationMultidim.h>` for
Powell, Nelder-Mead, and quasi-Newton methods, and
`<mml/algorithms/Optimization/LinearProgramming.h>` for simplex-based LP models.

### 5.1 Real function — golden section and Brent

**Problem:** minimize a smooth real function of one variable without supplying derivatives.
Golden section is simple and robust; Brent usually reaches the same minimum in fewer iterations
when the function is smooth enough for parabolic interpolation.

**Setup**

```cpp
#include <mml/algorithms/Optimization/Optimization.h>

RealFunction objective([](Real x) {
    return (x - 1.25) * (x - 1.25) + 0.15 * std::sin(3.0 * x);
});
```

**Solve**

```cpp
MinimumBracket bracket = Minimization::BracketMinimum(objective, 0.0, 1.0);

MinimizationResult golden = Minimization::GoldenSectionSearch(objective, bracket, 1e-10);
MinimizationResult brent = Minimization::BrentMinimize(objective, bracket, 1e-10);

// golden.xmin and brent.xmin should agree; brent.iterations is usually smaller here.
```

Bracket first when you want to compare methods fairly. Both methods then solve the same local
minimization problem inside the same bracket.

### 5.2 Multidimensional function — Powell, Nelder-Mead, and BFGS

**Problem:** minimize a two-variable objective from a starting guess. Powell and Nelder-Mead are
derivative-free; BFGS uses gradient information when it is available.

**Setup**

```cpp
#include <mml/algorithms/Optimization/OptimizationMultidim.h>

class OffsetBowl : public Optimization::IDifferentiableScalarFunction<2> {
public:
    Real operator()(const VectorN<Real, 2>& x) const override
    {
        Real dx = x[0] - 1.0;
        Real dy = x[1] + 2.0;
        return dx * dx + 2.0 * dy * dy;
    }

    void Gradient(const VectorN<Real, 2>& x, VectorN<Real, 2>& grad) const override
    {
        grad[0] = 2.0 * (x[0] - 1.0);
        grad[1] = 4.0 * (x[1] + 2.0);
    }
};

OffsetBowl objective;
VectorN<Real, 2> start{ -3.0, 4.0 };
```

**Solve**

```cpp
auto powell = Optimization::PowellMinimize<2>(objective, start);
auto nelderMead = Optimization::NelderMeadMinimize<2>(objective, start, 0.75, 1e-10);
auto bfgs = Optimization::BFGSMinimize<2>(objective, start);

// All three should converge near x = (1, -2), f = 0.
```

Use Powell as a strong derivative-free default for smooth moderate-dimensional objectives.
Nelder-Mead is convenient when only function values are trustworthy. Use BFGS when gradients are
available and reasonably accurate.

### 5.3 Linear programming example

**Problem:** choose production amounts for two products to maximize profit subject to resource
constraints. Variables are nonnegative by default.

**Setup**

```cpp
#include <mml/algorithms/Optimization/LinearProgramming.h>

Optimization::LinearProgram lp(2, "ProductionPlan");
lp.SetVariableNames({ "standard", "premium" });
lp.SetObjective({ 40.0, 30.0 }, Optimization::LPObjective::Maximize);

lp.AddConstraint({ 2.0, 1.0 }, Optimization::LPConstraintType::LessEqual,
                 100.0, "machine-hours");
lp.AddConstraint({ 1.0, 1.0 }, Optimization::LPConstraintType::LessEqual,
                 80.0, "assembly-hours");
lp.AddConstraint({ 1.0, 0.0 }, Optimization::LPConstraintType::LessEqual,
                 40.0, "standard-demand");
```

**Solve**

```cpp
Optimization::LPResult result = Optimization::SolveLP(lp);

if (result.IsOptimal()) {
    Real standard = result.x[0];
    Real premium = result.x[1];
    Real profit = result.objectiveValue;
}
```

For this model the optimum is 20 standard units and 60 premium units, with profit 2600. Slack
values, dual values, and reduced costs are available on `LPResult` for diagnostics.

**Notes**

- Bracketing matters for one-dimensional minimization; it identifies the local minimum that
  golden section and Brent will refine.
- Multidimensional methods return `xmin`, `fmin`, iteration counts, and convergence status.
- Dense LP examples are small by design; large production LPs usually need sparse/revised
  simplex tooling.
- Runnable version: `Cookbook_Recipe05_Optimization()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 6: Solving differential equations

Systems implement `IODESystem` (`getDim()` + `derivs(t, y, dydt)`); the quickest way is the
`ODESystem` adapter over a lambda. Solvers live in `<mml/algorithms/ODESolvers.h>`
(+ `ODESolvers/ODESolverEventDetection.h` for events).

### 6.1 Fixed-step RK4 on a simple problem

**Problem:** harmonic oscillator x'' = −x over one period — the solution must return
exactly to its starting point (1, 0).

```cpp
ODESystem sho(2, [](Real t, const Vector<Real>& y, Vector<Real>& dydt) {
    dydt[0] = y[1];
    dydt[1] = -y[0];
});
Vector<Real> y0({ 1.0, 0.0 });

RungeKutta4_StepCalculator rk4;
ODESystemFixedStepSolver solver(sho, rk4);
auto sol = solver.integrate(y0, 0.0, 2 * Constants::PI, 100);   // 100 fixed steps
// sol.getXValue(100, 0) = 0.99999996 - seven digits from 100 RK4 steps
```

Other step calculators: `EulerStep_Calculator`, `VelocityVerlet_StepCalculator`
(Hamiltonian systems), `RK5_CashKarp_StepCalculator`.

### 6.2 Adaptive Dormand-Prince on a hard problem

**Problem:** the chaotic Lorenz attractor — trajectories bend violently, so the step size
must adapt.

```cpp
ODESystem lorenz(3, [](Real t, const Vector<Real>& y, Vector<Real>& dydt) {
    dydt[0] = 10.0 * (y[1] - y[0]);
    dydt[1] = y[0] * (28.0 - y[2]) - y[1];
    dydt[2] = y[0] * y[1] - 8.0 / 3.0 * y[2];
});
Vector<Real> y0({ 1.0, 1.0, 1.0 });

ODEAdaptiveIntegrator<DormandPrince5_Stepper> integrator(lorenz);
auto sol = integrator.integrate(y0, 0.0, 10.0, 0.1, 1e-8);   // save every 0.1, eps = 1e-8
// sol.size() saved points; dense output via sol.getSolAsSplineInterp(comp)
```

Steppers: `CashKarp_Stepper`, `DormandPrince5_Stepper` (FSAL), `DormandPrince8_Stepper`
(high accuracy); aliases `CashKarpIntegrator`, `DormandPrince5Integrator`, ...

### 6.3 Event detection — artillery shot

**Problem:** a shell fired at 45°, 100 m/s, with quadratic air drag — integrate until it
hits the ground. “When does y cross zero going down?” is an **event**, not a fixed end time.

```cpp
#include <mml/algorithms/ODESolvers/ODESolverEventDetection.h>

class ArtilleryShot : public IODESystemWithEvents {
    Real _g = 9.81, _drag = 0.001;
public:
    int getDim() const override { return 4; }   // (x, y, vx, vy)

    void derivs(Real t, const Vector<Real>& s, Vector<Real>& dsdt) const override {
        Real speed = std::sqrt(s[2] * s[2] + s[3] * s[3]);
        dsdt[0] = s[2];
        dsdt[1] = s[3];
        dsdt[2] = -_drag * speed * s[2];          // quadratic air drag
        dsdt[3] = -_g - _drag * speed * s[3];
    }

    int getNumEvents() const override { return 1; }
    Real eventFunction(int, Real, const Vector<Real>& s) const override { return s[1]; } // y = 0
    EventDirection getEventDirection(int) const override { return EventDirection::Decreasing; }
    EventAction getEventAction(int) const override { return EventAction::Stop; }
};

ArtilleryShot shot;
Vector<Real> s0({ 0.0, 0.0, 70.71, 70.71 });   // 45 degrees at 100 m/s

DormandPrince5EventIntegrator integrator(shot);
auto result = integrator.integrateWithEvents(shot, s0, 0.0, 60.0, 0.05, 1e-10, 1e-12);

// result.terminatedByEvent == true; the impact event pinpoints landing:
// t = 12.251 s, range = 596.1 m   (the vacuum range would be v0^2/g = 1019.4 m)
const auto& impact = result.events.back();     // impact.time, impact.state
```

Events can also `Restart` with a modified state (bouncing ball) or just `Continue`
(recording apex crossings), with direction filtering per event.

### 6.4 Boundary value problem — shooting method

**Problem:** y'' = −y with **boundary** conditions y(0) = 0 and y(π/2) = 1 (not initial
ones!). Exact solution y = sin t, so the missing initial slope is y'(0) = 1. The shooting
method finds it by root-finding on the boundary mismatch:

```cpp
ODESystem sho(2, [](Real t, const Vector<Real>& y, Vector<Real>& dydt) {
    dydt[0] = y[1];
    dydt[1] = -y[0];
});

BVPShootingSolver solver(sho);
auto result = solver.solvePositionBVP(0.0, Constants::PI / 2,
                                      0.0, 1.0,     // y(0) = 0, y(pi/2) = 1
                                      0.5, 2.0);    // two guesses for y'(0)

// result.converged, result.shootingParams[0] = 1.0 (machine precision),
// result.residual ~ 1e-16, result.solution = full trajectory
```

`solve1D` handles general component/target combinations; `BVPShootingSolverND` shoots
multiple parameters at once (e.g., 3D projectile aiming).

### 6.5 Where to aim? — shooting with a free boundary

**Problem:** the gun from §6.3 fires at a fixed muzzle speed v₀ = 100 m/s. **At what
elevation must it shoot to land exactly 400 m away?** This is the *original* shooting
problem — and unlike §6.4 the “boundary” (the landing) happens at an unknown time, so the
fixed-interval BVP solver doesn't apply. Instead, compose two recipes: the impact *event*
(§6.3) gives range as a function of elevation, and *root finding* (Recipe 2) aims the gun:

```cpp
const Real v0 = 100.0, targetRange = 400.0;
const Real deg = Constants::PI / 180.0;

// range(theta): fire the ArtilleryShot system from 6.3, read the impact event
auto rangeFor = [&](Real theta) {
    ArtilleryShot shot;
    Vector<Real> s0({ 0.0, 0.0, v0 * std::cos(theta), v0 * std::sin(theta) });
    DormandPrince5EventIntegrator integrator(shot);
    auto res = integrator.integrateWithEvents(shot, s0, 0.0, 60.0, 0.05, 1e-10, 1e-12);
    return res.events.back().state[0];
};

// Root-find the elevation on the miss distance. With drag there are TWO solutions:
auto miss = [&](Real theta) { return rangeFor(theta) - targetRange; };

RootFinding::RootFindingConfig config;
config.tolerance = 1e-10;
auto lowSol  = RootFinding::FindRootBrent(miss,  5 * deg, 45 * deg, config);  // direct fire
auto highSol = RootFinding::FindRootBrent(miss, 45 * deg, 85 * deg, config);  // mortar lob

// direct fire at 15.82 deg, mortar lob at 68.13 deg - both land at 400.000 m
```

Every Brent iteration fires a full drag-affected trajectory and locates its impact by
event bisection — numerical artillery aiming in ~25 lines, with no closed-form solution
anywhere in sight.

### 6.6 Stiff system — Rosenbrock 2(3)

**Problem:** y' = −1000 (y − sin t) + cos t, y(0) = 0 — exact solution y = sin t, but the
stiffness constant 1000 would force an explicit method to h < 0.002 for *stability* even
though the solution is smooth. Stiff solvers need the Jacobian ∂f/∂y, so the system
implements `IODESystemWithJacobian`:

```cpp
class StiffRelaxation : public IODESystemWithJacobian {
public:
    int getDim() const override { return 1; }

    void derivs(Real t, const Vector<Real>& y, Vector<Real>& dydt) const override {
        dydt[0] = -1000.0 * (y[0] - std::sin(t)) + std::cos(t);
    }

    void jacobian(const Real t, const Vector<Real>& y, Vector<Real>& dydt, Matrix<Real>& J) const override {
        J(0, 0) = -1000.0;
    }
};

StiffRelaxation stiff;
Vector<Real> y0({ 0.0 });

Rosenbrock23Solver solver(stiff, 1e-6, 1e-6);   // 2nd-order: keep tolerances moderate
auto sol = solver.Solve(0.0, y0, 1.0, 0.1);
// y(1) = 0.8414694 vs exact sin 1 = 0.8414710
```

Siblings for stiff problems: `SolveBackwardEuler`, BDF2, and the DAE solver family
(`mml/algorithms/DAESolvers/`) for the truly hard cases.

**Notes**

- Fixed-step integrators return the trajectory on your grid; adaptive ones control local
  error and can interpolate (`getSolAsSplineInterp`) for dense output.
- Event times are located by bisection to the event tolerance — far more accurate than
  reading off the nearest saved point.
- Runnable version: `Cookbook_Recipe06_DifferentialEquations()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 7: Eigenvalues and eigenvectors

All four cases use the aggregate header `<mml/algorithms/EigenSystemSolvers.h>`. Choose the
solver from the matrix structure: symmetric and Hermitian solvers exploit their stronger
mathematical guarantees, while the general solvers allow complex spectra.

### 7.1 Real symmetric matrix

**Problem:** find all eigenvalues and eigenvectors of a real symmetric matrix. Symmetry
guarantees real eigenvalues and an orthonormal basis of eigenvectors.

**Setup** — this matrix has eigenvalues 3, 3, and 6:

```cpp
#include <mml/algorithms/EigenSystemSolvers.h>

Matrix<Real> A(3, 3, { 4.0, 1.0, 1.0,
                       1.0, 4.0, 1.0,
                       1.0, 1.0, 4.0 });
```

**Solve** — eigenvalues are sorted ascending and column `i` of `eigenvectors` belongs to
`eigenvalues[i]`:

```cpp
auto result = SymmMatEigenSolverJacobi::Solve(A);

for (int i = 0; i < A.cols(); i++) {
    Vector<Real> v = result.eigenvectors.VectorFromColumn(i);
    Vector<Real> residual = A * v - result.eigenvalues[i] * v;
    // residual.NormL2() is near machine precision
}
// result.eigenvalues = (3, 3, 6), result.converged == true
```

For large real symmetric matrices, `SymmMatEigenSolverQR` has the same result shape and is
usually faster; Jacobi is a robust choice for small and medium dense matrices.

### 7.2 General real matrix

**Problem:** find the eigenvalues of a nonsymmetric real matrix. Unlike the symmetric case,
the eigenvalues need not be real. A planar 90-degree rotation supplies the conjugate pair
$+i,-i$, while the third axis is scaled by 2.

**Setup**

```cpp
Matrix<Real> A(3, 3, { 0.0, -1.0, 0.0,
                       1.0,  0.0, 0.0,
                       0.0,  0.0, 2.0 });
```

**Solve**

```cpp
auto result = EigenSolver::Solve(A);

for (const auto& lambda : result.eigenvalues)
    std::cout << lambda << "\n";
// eigenvalues: +i, -i, 2 (order is not guaranteed)
```

Each value has `.real`, `.imag`, `.isComplex()`, and `.magnitude()`. The solver also returns
eigenvectors, encoding a conjugate pair's real and imaginary parts in adjacent real columns;
when only the spectrum is needed, the `eigenvalues` collection is sufficient.

### 7.3 Complex Hermitian matrix

**Problem:** find all eigenvalues and eigenvectors of a complex Hermitian matrix $A=A^*$.
Hermitian matrices are the complex analogue of real symmetric matrices: all eigenvalues are
real and the eigenvectors can be chosen unitary.

**Setup**

```cpp
Matrix<Complex> A(3, 3, {
    {2.0, 0.0},  {1.0, 1.0},  {0.0, -1.0},
    {1.0, -1.0}, {3.0, 0.0},  {2.0, 0.5},
    {0.0, 1.0},  {2.0, -0.5}, {5.0, 0.0}
});
```

**Solve**

```cpp
auto result = HermitianMatEigenSolverJacobi::Solve(A);

for (int i = 0; i < A.cols(); i++) {
    Vector<Complex> v = result.eigenvectors.VectorFromColumn(i);
    Vector<Complex> residual = A * v - Complex(result.eigenvalues[i], 0.0) * v;
    // residual.NormL2() is near machine precision
}
// result.eigenvalues is Vector<Real>, sorted ascending
```

The solver validates the Hermitian contract, including real diagonal entries, rather than
silently replacing the input with its Hermitian part.

### 7.4 General complex matrix

**Problem:** find the eigenvalues of a general complex matrix with no symmetry guarantee.
For a triangular matrix the answer is visible on the diagonal, making a useful check of the
general solver.

**Setup**

```cpp
Matrix<Complex> A(3, 3, {
    {1.0, 1.0}, {2.0, 0.0},  {0.0, 0.0},
    {0.0, 0.0}, {-2.0, 0.5}, {1.0, -1.0},
    {0.0, 0.0}, {0.0, 0.0},  {3.0, -2.0}
});
```

**Solve**

```cpp
auto result = ComplexEigenSolver::Solve(A);

for (Complex lambda : result.eigenvalues)
    std::cout << lambda << "\n";
// eigenvalues: 1+i, -2+0.5i, 3-2i
// result.maxResidual measures the worst returned right eigenpair
```

`ComplexEigenSolver` uses shifted complex QR and returns complex right eigenvectors as
columns. Check `result.converged` before consuming the decomposition and use
`result.maxResidual` as the quick numerical quality diagnostic.

**Notes**

- Prefer the structured solver whenever symmetry or Hermitian structure is known; it gives
  stronger guarantees and avoids unnecessary general-complex arithmetic.
- Eigenvectors are stored in columns in all four result types.
- Repeated eigenvalues can have different but equally valid eigenvector bases, so verify
  $Av=\lambda v$ instead of comparing vector components directly.
- Runnable version: `Cookbook_Recipe07_EigenvaluesAndEigenvectors()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 8: Matrix decompositions

The stateless decomposition API lives in `<mml/algorithms/MatrixAlg.h>`. SVD is the most
diagnostic decomposition: it exposes rank, null spaces, conditioning, and controlled low-rank
approximations. QR, LU, and Cholesky are narrower tools for orthogonal bases, pivoted square
factorization, and symmetric positive-definite matrices.

### 8.1 SVD: singular values, rank, and condition number

**Problem:** inspect a rectangular data matrix and verify that its SVD really reconstructs the
input.

**Setup**

```cpp
#include <mml/algorithms/MatrixAlg.h>

Matrix<Real> A(4, 3, { 1.0,  0.0,  2.0,
                       0.0,  1.0, -1.0,
                       2.0,  1.0,  0.0,
                       1.0, -1.0,  1.0 });
```

**Solve**

```cpp
auto svd = MatrixAlg::SVDDecompose(A);

Matrix<Real> sigma(A.rows(), A.cols());
for (int i = 0; i < svd.singularValues.size(); ++i)
    sigma(i, i) = svd.singularValues[i];

Matrix<Real> reconstructed = svd.U * sigma * svd.V.transpose();
Real residual = MatrixAlg::FrobeniusNorm(A - reconstructed);
// svd.rank is the numerical rank.
// MatrixAlg::ConditionNumber(A) uses the same SVD machinery.
```

Use the residual as the quick sanity check. For a well-behaved dense matrix it should be close
to floating-point roundoff relative to `MatrixAlg::FrobeniusNorm(A)`.

### 8.2 SVD: low-rank approximation

**Problem:** compress a measurement matrix by keeping only the strongest singular directions.

**Setup**

```cpp
Matrix<Real> measurements(5, 4, { 10.0,  9.8,  6.0,  5.9,
                                  8.0,  7.9,  4.8,  4.7,
                                  6.0,  6.1,  3.5,  3.6,
                                  4.0,  4.1,  2.5,  2.4,
                                  2.0,  2.1,  1.2,  1.1 });
```

**Solve**

```cpp
auto svd = MatrixAlg::SVDDecompose(measurements);

Matrix<Real> sigma1(measurements.rows(), measurements.cols());
sigma1(0, 0) = svd.singularValues[0];
Matrix<Real> rank1 = svd.U * sigma1 * svd.V.transpose();

Matrix<Real> sigma2 = sigma1;
sigma2(1, 1) = svd.singularValues[1];
Matrix<Real> rank2 = svd.U * sigma2 * svd.V.transpose();

Real rank1Error = MatrixAlg::FrobeniusNorm(measurements - rank1);
Real rank2Error = MatrixAlg::FrobeniusNorm(measurements - rank2);
```

Adding a singular value cannot increase the best Frobenius-norm approximation error; the second
error should be smaller or equal.

### 8.3 SVD: null space and rank deficiency

**Problem:** find directions `x` where `A*x = 0` for a matrix with dependent rows.

**Setup**

```cpp
Matrix<Real> A(3, 4, { 1.0, 2.0, 0.0,  3.0,
                       2.0, 4.0, 0.0,  6.0,
                       0.0, 1.0, 1.0, -1.0 });
```

**Solve**

```cpp
auto spaces = MatrixAlg::FundamentalSubspacesOf(A);
Matrix<Real> N = spaces.nullSpace;          // columns form an orthonormal basis for null(A)
Real nullResidual = MatrixAlg::FrobeniusNorm(A * N);
// spaces.rank == 2, N.cols() == 2, nullResidual is near zero
```

`FundamentalSubspacesOf` also returns column space, row space, and left null space bases from
the same SVD pass.

### 8.4 QR: orthogonal basis and reconstruction

**Problem:** turn the columns of a tall matrix into an orthonormal basis and triangular factor.

**Setup**

```cpp
Matrix<Real> A(4, 2, { 1.0, 2.0,
                       2.0, 0.0,
                       0.0, 3.0,
                       1.0, 1.0 });
```

**Solve**

```cpp
auto qr = MatrixAlg::QRDecompose(A);
Real reconstruction = MatrixAlg::FrobeniusNorm(A - qr.Q * qr.R);
Real orthogonality = MatrixAlg::FrobeniusNorm(
    qr.Q.transpose() * qr.Q - Matrix<Real>::Identity(qr.Q.cols()));
```

The MML QR helper returns economy QR: `Q` has one column per input column, and `R` is square.

### 8.5 LU: factors and determinant

**Problem:** factor a square matrix with pivoting, then verify the row-permuted factorization.

**Setup**

```cpp
Matrix<Real> A(3, 3, { 0.0, 2.0, 1.0,
                       1.0, 1.0, 0.0,
                       2.0, 0.0, 1.0 });
```

**Solve**

```cpp
auto lu = MatrixAlg::LUDecompose(A);

Matrix<Real> PA(A.rows(), A.cols());
for (int row = 0; row < A.rows(); ++row)
    for (int col = 0; col < A.cols(); ++col)
        PA(row, col) = A(lu.permutation[row], col);

Real residual = MatrixAlg::FrobeniusNorm(lu.L * lu.U - PA);
// lu.determinant is the determinant with pivot signs accounted for
```

Use LU when you want a square factorization that can be reused for determinants, inverses, and
linear solves.

### 8.6 Cholesky: positive-definite decomposition

**Problem:** factor a symmetric positive-definite covariance-like matrix as $A=L L^T$.

**Setup**

```cpp
Matrix<Real> covariance(3, 3, { 4.0, 2.0, 0.6,
                                2.0, 3.0, 0.5,
                                0.6, 0.5, 1.5 });
```

**Solve**

```cpp
auto chol = MatrixAlg::CholeskyDecompose(covariance);
Real residual = MatrixAlg::FrobeniusNorm(chol.L * chol.L.transpose() - covariance);
// chol.L is lower triangular
```

Cholesky is faster and more structure-preserving than LU for positive-definite matrices, but
it should fail on matrices that do not satisfy that contract.

**Notes**

- SVD is the safest exploratory decomposition when rank or conditioning is uncertain.
- QR requires rows >= columns in the economy helper.
- LU and Cholesky are square-matrix tools; Cholesky additionally requires positive definiteness.
- Runnable version: `Cookbook_Recipe08_MatrixDecompositions()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 9: Matrix properties and diagnostics

The property and diagnostic helpers are also in `<mml/algorithms/MatrixAlg.h>`. Use the
stateless functions for one-off checks, or `MatrixAnalyzer` from
`<mml/algorithms/Analyzers/MatrixAnalyzer.h>` when several related properties should be cached
against the same matrix.

### 9.1 Inspect matrix structure

**Problem:** classify a matrix before choosing a solver or decomposition.

**Setup**

```cpp
#include <mml/algorithms/Analyzers/MatrixAnalyzer.h>

Matrix<Real> A(3, 3, { 4.0, 1.0, 0.0,
                       1.0, 3.0, 0.0,
                       0.0, 0.0, 2.0 });
```

**Inspect**

```cpp
MatrixAnalyzer<Real> analyzer(A);
bool square = analyzer.IsSquare();
bool symmetric = analyzer.IsSymmetric();
bool diagonal = analyzer.IsDiagonal();
bool dominant = analyzer.IsDiagonallyDominant();
```

This kind of check is useful before deciding between generic, symmetric, and positive-definite
algorithms.

### 9.2 Norms and approximation error

**Problem:** measure matrix scale and quantify the error of an approximate result.

**Setup**

```cpp
Matrix<Real> A(3, 3, { 3.0, -1.0, 0.0,
                       2.0,  4.0, 1.0,
                       0.0, -2.0, 5.0 });
Matrix<Real> approximation(3, 3, { 3.0, -1.0, 0.0,
                                   2.0,  4.0, 0.9,
                                   0.0, -2.0, 5.0 });
```

**Measure**

```cpp
Real frobenius = MatrixAlg::FrobeniusNorm(A);
Real oneNorm = MatrixAlg::OneNorm(A);
Real infinityNorm = MatrixAlg::InfinityNorm(A);
Real error = MatrixAlg::FrobeniusNorm(A - approximation);
```

The Frobenius norm is often the most convenient dense-matrix error measure; the one and infinity
norms expose maximum column and row sums.

### 9.3 Trace, determinant, and rank

**Problem:** compute compact invariants for a square matrix.

**Setup**

```cpp
Matrix<Real> A(3, 3, { 2.0, 1.0, 0.0,
                       1.0, 3.0, 1.0,
                       0.0, 1.0, 2.0 });
```

**Measure**

```cpp
Real trace = MatrixAlg::Trace(A);
Real determinant = MatrixAlg::Determinant(A);
int rank = MatrixAlg::Rank(A);              // SVD-based
int gaussianRank = MatrixAlg::RankGaussian(A);
```

Use SVD rank when numerical rank matters; Gaussian rank is a direct algebraic diagnostic with a
fixed tolerance.

### 9.4 Condition number and sensitivity

**Problem:** decide whether a matrix is likely to amplify input or roundoff errors.

**Setup**

```cpp
Matrix<Real> A(2, 2, { 1.0, 0.999,
                       0.999, 0.998 });
```

**Assess**

```cpp
MatrixAnalyzer<Real> analyzer(A);
Real cond2 = analyzer.ConditionNumber();
auto stability = analyzer.AssessStability();
auto digitsLost = analyzer.ExpectedDigitsLost();
```

Large condition numbers do not make every result wrong, but they tell you how much care to take
with residual checks and data precision.

### 9.5 Inverse and identity verification

**Problem:** compute an inverse only when you really need the inverse matrix, then verify it.

**Setup**

```cpp
Matrix<Real> A(3, 3, { 4.0, 1.0, 0.0,
                       1.0, 3.0, 1.0,
                       0.0, 1.0, 2.0 });
```

**Verify**

```cpp
Matrix<Real> inverse = MatrixAlg::Inverse(A);
Real residual = MatrixAlg::FrobeniusNorm(A * inverse - Matrix<Real>::Identity(A.rows()));
```

For solving `A*x=b`, prefer a linear solver or decomposition reuse. The inverse is best reserved
for diagnostics, explicit transformations, and formulas that genuinely require it.

### 9.6 Matrix slicing and transformations

**Problem:** extract a rectangular block and build a transformed matrix from it.

**Setup**

```cpp
Matrix<Real> samples(4, 4, { 1.0,  2.0,  3.0,  4.0,
                             5.0,  6.0,  7.0,  8.0,
                             9.0, 10.0, 11.0, 12.0,
                            13.0, 14.0, 15.0, 16.0 });
```

**Transform**

```cpp
Matrix<Real> center(samples, 1, 1, 2, 2);   // rows 1..2, cols 1..2
Matrix<Real> gram = center.transpose() * center;
Real centerTrace = MatrixAlg::Trace(center);
```

The submatrix constructor copies the requested block, so `center` can outlive the source matrix.

### 9.7 Matrix-valued functions

**Problem:** represent a matrix that depends on a scalar parameter and evaluate diagnostics at
particular parameter values.

**Setup**

```cpp
class TimeDependentMatrix : public IFunction<Matrix<Real>, Real>
{
public:
    Matrix<Real> operator()(Real t) const override
    {
        return Matrix<Real>(2, 2, { 1.0 + t, 0.25,
                                   0.25, 2.0 - 0.5 * t });
    }
};
```

**Evaluate**

```cpp
TimeDependentMatrix A;
Matrix<Real> atQuarter = A(0.25);
Matrix<Real> atHalf = A(0.5);
Matrix<Real> rate = (atHalf - atQuarter) * (1.0 / 0.25);
```

**Measure**

```cpp
Real trace = MatrixAlg::Trace(atHalf);
Real determinant = MatrixAlg::Determinant(atHalf);
```

The generic `IFunction<Return, Argument>` interface is enough for matrix-valued maps. Once a
matrix value is produced, all ordinary `MatrixAlg` diagnostics apply.

**Notes**

- Structure predicates help choose specialized algorithms before doing expensive work.
- Norms and residuals are the basic numerical sanity checks for approximations.
- `MatrixAnalyzer` is most useful when multiple diagnostics share a decomposition or other
  cached intermediate.
- Runnable version: `Cookbook_Recipe09_MatrixPropertiesAndDiagnostics()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 10: Statistics

The basic API is in `<mml/algorithms/Statistics.h>`. Histograms and probability distributions
use `<mml/algorithms/Statistics/Histogram.h>` and `<mml/algorithms/Statistics/Distributions.h>`.

### 10.1 Quick descriptive summary

**Problem:** summarize the center and spread of a dataset while retaining the distinction
between a sample and a complete population.

```cpp
Vector<Real> data({ 2, 4, 4, 4, 5, 5, 7, 9 });

Real mean, sampleStdDev, q1, median, q3, minimum, maximum;
Statistics::AvgStdDev(data, mean, sampleStdDev);
Statistics::Quartiles(data, q1, median, q3);
Statistics::MinMax(data, minimum, maximum);
Real populationStdDev = Statistics::PopulationStdDev(data);
// mean=5, median=4.5, sampleStdDev=2.13809, populationStdDev=2
```

`StdDev` and `SampleStdDev` divide variance by $n-1$; `PopulationStdDev` divides by $n$.

### 10.2 Robust statistics with outliers

**Problem:** one 100 ms timeout contaminates otherwise typical 10–14 ms response times.

```cpp
Vector<Real> responseMs({ 10, 12, 11, 13, 12, 11, 14, 100 });

Real mean        = Statistics::Mean(responseMs);              // 23.875
Real median      = Statistics::Median(responseMs);            // 12
Real stddev      = Statistics::StdDev(responseMs);            // 31.188
Real trimmedMean = Statistics::TrimmedMean(responseMs, 12.5); // 12.167
Real mad         = Statistics::MAD(responseMs);               // 1
Real iqr         = Statistics::IQR(responseMs);               // 2.25
```

The mean and standard deviation answer “including the timeout”; median, trimmed mean, MAD,
and IQR describe normal behavior without letting one extreme observation dominate.

### 10.3 Weighted measurements

**Problem:** calculate a course grade when subjects carry different credit weights, and measure
the weighted relationship between study time and grade.

```cpp
Vector<Real> grades({ 85, 90, 78, 92 });
Vector<Real> credits({ 2, 3, 1, 4 });
Vector<Real> studyHours({ 6, 8, 4, 10 });

Real grade = Statistics::WeightedMean(grades, credits); // 88.6
Real spread = Statistics::WeightedStdDev(grades, credits);
Real correlation = Statistics::WeightedPearsonCorrelation(studyHours, grades, credits);
```

Weights are frequency weights: weighted variance uses denominator $\sum w_i-1$.

### 10.4 Relationships between variables

**Problem:** quantify the relationship between study hours and test scores, then inspect several
variables together. Rows are observations; columns are variables.

```cpp
Vector<Real> studyHours({ 2, 3, 4, 5, 6, 7, 8, 9 });
Vector<Real> testScores({ 65, 70, 72, 80, 85, 88, 92, 95 });

Real covariance = Statistics::Covariance(studyHours, testScores);
Real r = Statistics::PearsonCorrelation(studyHours, testScores); // 0.9924
Real rSquared = Statistics::RSquared(studyHours, testScores);    // 0.9849

Matrix<Real> observations(8, 3);
for (int row = 0; row < 8; ++row) {
    observations(row, 0) = studyHours[row];
    observations(row, 1) = testScores[row];
    observations(row, 2) = 10.0 - 0.5 * studyHours[row];
}
Matrix<Real> correlations = Statistics::CorrelationMatrix(observations);
```

Correlation measures linear association, not causation. $R^2=r^2$ is the fraction of variance
explained by a one-predictor linear relationship.

### 10.5 Histogram and empirical CDF

**Problem:** understand delivery-time shape and answer “what fraction arrives within 30 minutes?”

```cpp
#include <mml/algorithms/Statistics/Histogram.h>

Vector<Real> deliveryMinutes({ 18, 19, 20, 21, 21, 22, 23, 24,
                               24, 25, 26, 27, 29, 31, 35, 48 });
auto histogram = Statistics::Histogram::ComputeHistogramAuto(
    deliveryMinutes, Statistics::Histogram::BinningMethod::FreedmanDiaconis);
auto cumulativeCounts = histogram.GetCumulativeCounts();
Real withinThirty = Statistics::Histogram::EvaluateECDF(deliveryMinutes, 30.0); // 0.8125
```

Freedman–Diaconis chooses bin width from sample size and IQR. The ECDF answers threshold
questions without assuming a theoretical distribution.

### 10.6 Confidence interval for an unknown mean

**Problem:** estimate mean package fill weight and its 95% confidence interval from ten
measurements when population variance is unknown.

```cpp
#include <mml/algorithms/Statistics/Distributions.h>

Vector<Real> fillWeights({ 500.2, 499.8, 500.5, 501.0, 499.4,
                           500.1, 500.7, 499.9, 500.3, 499.6 });
Real mean, sampleStdDev;
Statistics::AvgStdDev(fillWeights, mean, sampleStdDev);
Statistics::TDistribution studentT(fillWeights.size() - 1);
Real critical = studentT.inverseCdf(0.975);
Real margin = critical * sampleStdDev / std::sqrt(static_cast<Real>(fillWeights.size()));
// 95% CI: [499.7943, 500.5057] g
```

This assumes independent observations from an approximately normal population. Student's $t$
is used because the standard deviation is estimated from the sample.

### 10.7 Real dataset — Palmer Penguins

**Question:** which measurements distinguish penguin species, and how strongly are body mass and
flipper length related?

The real 344-observation dataset is loaded from
[`test_data/statistics/palmer_penguins.csv`](../test_data/statistics/palmer_penguins.csv).
Missing measurements remain in the source and are removed pairwise.

```cpp
#include <mml/tools/DataLoader.h>

Data::Dataset penguins = Data::LoadCSV("test_data/statistics/palmer_penguins.csv");
auto species = penguins.GetStringColumn("species");
const auto& massColumn = penguins["body_mass_g"];
const auto& flipperColumn = penguins["flipper_length_mm"];

Vector<Real> bodyMass, flipperLength;
std::map<std::string, Vector<Real>> massBySpecies;
for (std::size_t row = 0; row < penguins.NumRows(); ++row) {
    Real mass = massColumn.GetReal(row);
    Real flipper = flipperColumn.GetReal(row);
    if (!std::isfinite(mass) || !std::isfinite(flipper)) continue;
    bodyMass.push_back(mass);
    flipperLength.push_back(flipper);
    massBySpecies[species[row]].push_back(mass);
}
```

**Results:** 342 complete pairs give mean mass 4201.75 g, median 4050 g, IQR 1200 g,
MAD 600 g, and correlation $r=0.8712$. Mean masses are approximately 3701 g (Adelie),
3733 g (Chinstrap), and 5076 g (Gentoo). The demo also computes an 11-bin histogram,
$P(\text{mass}\le5000\text{ g})=0.8216$, and a 95% mean interval $[4116.46,4287.05]$ g.

The data were collected by Kristen Gorman and Palmer Station LTER and published by Horst,
Hill, and Gorman (2020), DOI [10.5281/zenodo.3960218](https://doi.org/10.5281/zenodo.3960218),
under CC0. Full provenance is in
[`test_data/statistics/README.md`](../test_data/statistics/README.md).

### 10.8 Real analysis — tyre degradation at Monza

**Question:** did tyre age measurably increase lap time during the 2024 Italian Grand Prix?

The checked-in [`lap_times.csv`](../src/examples/03_formula_1_sim/data/2024/16_italian_grand_prix/lap_times.csv)
contains driver, race lap, lap time, sector times, compound, and tyre life. The demo loads it with
`Data::LoadCSV`, detects stints from compound changes and tyre-life resets, rejects incomplete timing,
separates compounds and drivers, and applies per-stint MAD filtering to remove pit/safety-car/extreme laps.

A naive correlation is invalid because tyre age and race lap rise together while fuel mass falls.
The demo removes each driver's compound-specific baseline and fits

$$
\Delta t=\beta_{\rm tyre}\,\Delta(\text{tyre age})
       +\beta_{\rm race}\,\Delta(\text{race lap}).
$$

**Results:** medium-tyre age and race lap are collinear in the available one-stint sequences, so
their separate effects are **not identifiable**; their raw combined trend is about $-0.021$ s/lap
over 229 cleaned observations. Strategy variation makes the hard-tyre model identifiable over 661
observations: estimated tyre degradation is **+0.054 s per tyre lap**, while the race-lap effect is
**−0.067 s per lap**, consistent with cars becoming faster as fuel burns off.

The raw medium trend appears to improve with age, but cannot distinguish fuel burn from tyre wear.
Weather, traffic, track evolution, setup, and driver management remain unmodeled confounders, so
these are observational associations rather than causal estimates.

**Notes**

- Validate missing values before passing columns to statistics functions.
- Use robust summaries for operational data containing timeouts, pit stops, or sensor faults.
- A confidence interval quantifies sampling uncertainty under model assumptions; it does not prove
  the assumptions or practical significance.
- Runnable version: `Cookbook_Recipe10_Statistics()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 11: 3D geometry: planes, lines, and bodies

The base 3D geometry primitives are in `<mml/base/Geometry/Geometry3D.h>`. Solid bodies and
bounding volumes are in `<mml/base/Geometry/Geometry3DBodies.h>`.

This recipe focuses on constructing and querying concrete geometry objects. Use Recipe 12 for
higher-level geometry algorithms built on these primitives.

### 11.1 Build points, vectors, lines, and segments

**Problem:** describe a surveyed baseline between two measured 3D stations, then evaluate the
corresponding infinite line and finite segment.

**Setup**

```cpp
Pnt3Cart stationA(0.0, 0.0, 0.0);
Pnt3Cart stationB(3.0, 4.0, 0.0);
Vec3Cart baseline(stationA, stationB);
Line3D surveyLine(stationA, baseline);
SegmentLine3D surveySegment(stationA, stationB);
```

**Use**

```cpp
Real length = baseline.NormL2();                    // 5
Vec3Cart direction = baseline.GetAsUnitVector();    // (0.6,0.8,0)
Real segmentLength = surveySegment.Length();        // 5
Pnt3Cart pointAfterTwoUnits = surveyLine(2.0);      // stationA + 2*baseline
```

`Line3D` represents the infinite parametric line. `SegmentLine3D` keeps the endpoints and clamps
distance queries to the finite segment.

### 11.2 Intersect a line with a plane

**Problem:** find where a descending path crosses the ground plane and detect a parallel miss.

**Setup**

```cpp
Plane3D ground = Plane3D::GetXYPlane();
Line3D descendingPath(Pnt3Cart(2.0, -1.0, 5.0), Vec3Cart(0.0, 1.0, -2.0));
Line3D parallelPath(Pnt3Cart(0.0, 0.0, 3.0), Vec3Cart(1.0, 0.0, 0.0));
```

**Solve**

```cpp
Pnt3Cart hitPoint;
bool hitsGround = ground.IntersectionWithLine(descendingPath, hitPoint);
bool parallelHitsGround = ground.IntersectionWithLine(parallelPath, hitPoint);
```

`IntersectionWithLine` returns `false` when the line is parallel to the plane. Prefer named plane
factories or point-and-normal construction over ambiguous coefficient-style construction.

### 11.3 Project points and measure distances

**Problem:** drop a sensor position onto an elevated platform and measure its distance from a
cable axis.

**Setup**

```cpp
Plane3D platform(Pnt3Cart(0.0, 0.0, 1.0), Vec3Cart(0.0, 0.0, 1.0));
Pnt3Cart sensor(2.0, 3.0, 5.0);
Line3D cableAxis(Pnt3Cart(0.0, 0.0, 0.0), Vec3Cart(1.0, 1.0, 0.0));
```

**Measure**

```cpp
Pnt3Cart footprint = platform.ProjectionToPlane(sensor);
Real heightAbovePlatform = platform.DistToPoint(sensor);

Pnt3Cart nearestCablePoint = cableAxis.NearestPointOnLine(sensor);
Real distanceToCable = cableAxis.Dist(sensor);
```

Plane projection uses the plane normal; line distance uses the perpendicular projection onto the
infinite line.

### 11.4 Work with triangles and rectangular surfaces

**Problem:** compute common properties of a triangular brace and a rectangular panel.

**Setup**

```cpp
Triangle3D brace(Pnt3Cart(0.0, 0.0, 0.0), Pnt3Cart(4.0, 0.0, 0.0), Pnt3Cart(0.0, 3.0, 0.0));

RectSurface3D panel(Pnt3Cart(0.0, 0.0, 0.0), Pnt3Cart(2.0, 0.0, 0.0),
                    Pnt3Cart(2.0, 3.0, 0.0), Pnt3Cart(0.0, 3.0, 0.0));
```

**Query**

```cpp
Real braceArea = brace.Area();
Pnt3Cart centroid = brace.Centroid();
bool loadPointInside = brace.IsPointInside(Pnt3Cart(1.0, 1.0, 0.0));

Real panelArea = panel.getArea();
VectorN<Real, 3> panelCenter = panel(0.0, 0.0);
```

`Triangle3D::IsPointInside` uses barycentric coordinates. `RectSurface3D` is a parametric surface
whose local coordinates are centered on the rectangle.

### 11.5 Create and query 3D bodies

**Problem:** model simple solid objects and query their volume, surface area, containment, and
centers.

**Setup**

```cpp
Cube3D equipmentBox(2.0, Pnt3Cart(1.0, 2.0, 3.0));
Sphere3D safetyZone(1.5, Pnt3Cart(1.0, 2.0, 3.0));
Cylinder3D tank(1.0, 3.0, Pnt3Cart(0.0, 0.0, 0.0));

Pnt3Cart probe(1.5, 2.0, 3.0);
```

**Query**

```cpp
Real cubeVolume = equipmentBox.Volume();
Real cubeArea = equipmentBox.SurfaceArea();
bool probeInBox = equipmentBox.IsInside(probe);

Real sphereVolume = safetyZone.Volume();
bool probeInZone = safetyZone.IsInside(probe);
Pnt3Cart tankCenter = tank.GetCenter();
```

Bodies share the `IBody` interface: `Volume`, `SurfaceArea`, `IsInside`, `GetCenter`,
`GetBoundingBox`, and `GetBoundingSphere`.

### 11.6 Use bounding volumes for spatial checks

**Problem:** use cheap bounding-volume tests before doing more expensive body-level work.

**Setup**

```cpp
Cube3D equipmentBox(2.0, Pnt3Cart(1.0, 2.0, 3.0));
Box3D bounds = equipmentBox.GetBoundingBox();
BoundingSphere3D sphereBounds = equipmentBox.GetBoundingSphere();
BoundingSphere3D nearbyObject(Pnt3Cart(3.0, 2.0, 3.0), 0.75);
```

**Check**

```cpp
Pnt3Cart probe(1.8, 2.0, 3.0);
Box3D expandedBounds = bounds.Expanded(0.5);

bool containsProbe = bounds.Contains(probe);
bool sphereOverlap = sphereBounds.Intersects(nearbyObject);
Real expandedVolume = expandedBounds.Volume();
```

Axis-aligned boxes are useful for fast containment and overlap checks. Bounding spheres are often
looser, but their overlap test is very cheap.

**Notes**

- Construct planes from a point and normal, three points, or the named helpers when possible.
- Use `Line3D` for infinite axes and `SegmentLine3D` when endpoint clamping matters.
- Bodies expose both exact body queries and conservative bounding volumes.
- Runnable version: `Cookbook_Recipe11_Geometry3D()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 12: 3D geometry: algorithms

The computational geometry aggregate header is `<mml/algorithms/ComputationalGeometry.h>`.
Point-cloud queries also use `<mml/algorithms/CompGeometry/KDTree.h>`.

This recipe covers the practical 3D algorithms currently present in MML: line classification,
closest approach, projection, plane intersections, triangle tests, ray picking, point-cloud
queries, and convex hulls. Full mesh processing, BVH acceleration, triangle-triangle intersection,
and 3D Delaunay tetrahedralization are deliberately outside this recipe for now.

### 12.1 Classify two 3D lines

**Problem:** determine whether two construction axes cross, are parallel, are coincident, or are
skew.

**Setup**

```cpp
Line3D craneRail(Pnt3Cart(0.0, 0.0, 0.0), Vec3Cart(1.0, 0.0, 0.0));
Line3D walkway(Pnt3Cart(2.0, -1.0, 0.0), Vec3Cart(0.0, 1.0, 0.0));
Line3D parallelRail(Pnt3Cart(0.0, 2.0, 0.0), Vec3Cart(1.0, 0.0, 0.0));
Line3D sameRail(Pnt3Cart(3.0, 0.0, 0.0), Vec3Cart(-1.0, 0.0, 0.0));
Line3D skewPipe(Pnt3Cart(0.0, 1.0, 1.0), Vec3Cart(0.0, 1.0, 0.0));
```

**Solve**

```cpp
auto crossing = craneRail.Intersection(walkway);
auto parallel = craneRail.Intersection(parallelRail);
auto coincident = craneRail.Intersection(sameRail);
auto skew = craneRail.Intersection(skewPipe);

Pnt3Cart crossingPoint = crossing.Point();
LineIntersectionType3D skewType = skew.type;
```

`Line3D::Intersection` returns a `LineIntersection3D` classification and, when meaningful, the
intersection point or closest-approach data.

### 12.2 Find closest approach between skew lines

**Problem:** compute the shortest separation between two 3D axes that do not intersect.

**Setup**

```cpp
Line3D cameraRay(Pnt3Cart(0.0, 0.0, 1.0), Vec3Cart(1.0, 0.0, 0.0));
Line3D robotAxis(Pnt3Cart(0.0, 2.0, 0.0), Vec3Cart(0.0, 0.0, 1.0));
```

**Solve**

```cpp
auto approach = cameraRay.Intersection(robotAxis);

Real clearance = approach.distance;
Pnt3Cart pointOnRay = approach.point1;
Pnt3Cart pointOnAxis = approach.point2;
```

For skew lines, `point1` and `point2` are the endpoints of the shortest connecting segment.

### 12.3 Project points onto lines, segments, and planes

**Problem:** measure where a sensor sits relative to an infinite cable, a finite boom, and an
elevated deck.

**Setup**

```cpp
Pnt3Cart sensor(2.0, 3.0, 5.0);
Line3D cable(Pnt3Cart(0.0, 0.0, 1.0), Vec3Cart(1.0, 1.0, 0.0));
SegmentLine3D boom(Pnt3Cart(0.0, 0.0, 0.0), Pnt3Cart(4.0, 0.0, 0.0));
Plane3D deck(Pnt3Cart(0.0, 0.0, 2.0), Vec3Cart(0.0, 0.0, 1.0));
```

**Measure**

```cpp
Pnt3Cart cableFoot = cable.NearestPointOnLine(sensor);
Real lineDistance = cable.Dist(sensor);
Real segmentDistance = boom.Dist(sensor);

Pnt3Cart deckFoot = deck.ProjectionToPlane(sensor);
Real planeDistance = deck.DistToPoint(sensor);
```

Use the segment type when endpoints matter. A line distance is perpendicular to an infinite line;
a segment distance is clamped to the finite endpoints.

### 12.4 Intersect lines with planes and planes with planes

**Problem:** find where a drill path hits the floor and where a floor/wall pair meet.

**Setup**

```cpp
Plane3D floor = Plane3D::GetXYPlane();
Plane3D wall(Pnt3Cart(2.0, 0.0, 0.0), Vec3Cart(1.0, 0.0, 0.0));
Line3D drillPath(Pnt3Cart(2.0, 1.0, 3.0), Vec3Cart(0.0, 0.0, -1.0));
```

**Solve**

```cpp
Pnt3Cart floorHit;
bool lineHitsFloor = floor.IntersectionWithLine(drillPath, floorHit);

Line3D floorWallEdge;
bool planesIntersect = floor.IntersectionWithPlane(wall, floorWallEdge);
```

Line-plane and plane-plane methods return `false` for parallel/no-solution cases. The plane-plane
result is a `Line3D` edge.

### 12.5 Work with triangle geometry and point containment

**Problem:** compute basic facet properties and determine whether a load point lies inside the
triangle.

**Setup**

```cpp
Triangle3D facet(Pnt3Cart(0.0, 0.0, 0.0),
                 Pnt3Cart(4.0, 0.0, 0.0),
                 Pnt3Cart(0.0, 3.0, 0.0));
Pnt3Cart loadPoint(1.0, 1.0, 0.0);
```

**Query**

```cpp
Real area = facet.Area();
Pnt3Cart centroid = facet.Centroid();
bool containsLoad = facet.IsPointInside(loadPoint);
Plane3D supportPlane = facet.getDefinedPlane();
```

Triangle containment expects the point to lie in the triangle plane. If your point is off-plane,
project it first or treat the off-plane distance separately.

### 12.6 Ray-pick a triangle

**Problem:** shoot a picking ray at a triangular facet and recover the hit point and barycentric
coordinates.

**Setup**

```cpp
Triangle3D facet(Pnt3Cart(0.0, 0.0, 0.0),
                 Pnt3Cart(4.0, 0.0, 0.0),
                 Pnt3Cart(0.0, 3.0, 0.0));
Pnt3Cart rayOrigin(1.0, 1.0, 5.0);
Vec3Cart rayDirection(0.0, 0.0, -1.0);
```

**Solve**

```cpp
auto hit = CompGeometry::Intersections::IntersectRayTriangle(rayOrigin, rayDirection, facet);
auto [u, v, w] = hit.GetBarycentricCoords();

bool wasHit = hit.hit;
Real rayParameter = hit.t;
Pnt3Cart hitPoint = hit.point;
```

The ray-triangle test is the standard building block for picking, visibility checks, and simple
collision probes. The barycentric coordinates are useful for interpolating per-vertex data at the
hit point.

### 12.7 Query a 3D point cloud efficiently

**Problem:** find nearby landmarks in a point cloud without scanning every point manually.

**Setup**

```cpp
std::vector<Pnt3Cart> landmarks = {
    Pnt3Cart(0.0, 0.0, 0.0), Pnt3Cart(2.0, 0.0, 0.0), Pnt3Cart(0.0, 2.0, 0.0),
    Pnt3Cart(0.0, 0.0, 2.0), Pnt3Cart(2.0, 2.0, 2.0), Pnt3Cart(5.0, 5.0, 0.0)
};
KDTree3D tree;
tree.build(landmarks);
```

**Query**

```cpp
Pnt3Cart query(1.1, 0.2, 0.1);

auto nearest = tree.findNearest(query);
auto neighbors = tree.findKNearest(query, 3);
auto inRadius = tree.findInRadius(query, 2.25);
auto inBox = tree.findInBox(Pnt3Cart(-0.1, -0.1, -0.1), Pnt3Cart(2.1, 2.1, 0.1));
```

Use nearest-neighbor queries for matching and snapping, radius queries for local neighborhoods, and
box queries for coarse spatial filtering.

### 12.8 Build a convex hull and test containment

**Problem:** construct the enclosing polyhedron of a 3D point set and use it as a volume envelope.

**Setup**

```cpp
std::vector<Pnt3Cart> points = {
    Pnt3Cart(0.0, 0.0, 0.0), Pnt3Cart(1.0, 0.0, 0.0),
    Pnt3Cart(1.0, 1.0, 0.0), Pnt3Cart(0.0, 1.0, 0.0),
    Pnt3Cart(0.0, 0.0, 1.0), Pnt3Cart(1.0, 0.0, 1.0),
    Pnt3Cart(1.0, 1.0, 1.0), Pnt3Cart(0.0, 1.0, 1.0),
    Pnt3Cart(0.4, 0.4, 0.4)
};
```

**Solve**

```cpp
auto hull = CompGeometry::ConvexHull3DComputer::Compute(points);

size_t vertices = hull.NumVertices();
size_t faces = hull.NumFaces();
Real volume = hull.Volume();
Real surfaceArea = hull.SurfaceArea();
bool containsCenter = hull.Contains(Pnt3Cart(0.5, 0.5, 0.5));
```

Interior points are ignored by the hull construction. Once built, the hull can be used for volume,
surface-area, centroid, face iteration, and containment queries.

**Notes**

- Use `Line3D::Intersection` when classification matters; use `Dist`/projection methods when you
  only need a distance.
- Prefer point-and-normal plane construction for readable code.
- `KDTree3D` is the practical tool for repeated point-cloud proximity queries.
- Convex hull containment is an envelope test, not a replacement for detailed mesh collision.
- Runnable version: `Cookbook_Recipe12_Geometry3DAlgorithms()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 13: Interpolating functions

The aggregate header is `<mml/base/InterpolatedFunction.h>`. It includes the 1D real-function
interpolants, 2D grid interpolants, and parametric curve interpolants used below.

Interpolation turns tabulated input data into callable functions. The practical choice is usually about
the shape of the data: linear interpolation is robust and local, low-order polynomial
interpolation is useful for smooth data and error estimates, cubic splines are the default for
smooth tabulated curves, and rational/barycentric forms help avoid high-degree polynomial
oscillation.

### 13.1 Interpolate a real function with several methods

**Problem:** you have a calibration table with measured values at non-uniform x positions. Build
interpolants from the two input vectors and estimate the value at an x that was not measured.

**Setup**

```cpp
#include <mml/base/InterpolatedFunction.h>

Vector<Real> x({ 0.0, 0.35, 0.8, 1.2, 1.65, 2.1, 2.55, 3.0, 3.45, 3.9, 4.4 });
Vector<Real> y({ 0.02, 0.28, 0.61, 0.86, 1.01, 0.92, 0.63, 0.21, -0.18, -0.52, -0.77 });

Real query = 2.35;
```

**Solve**

```cpp
LinearInterpRealFunc linear(x, y);
PolynomInterpRealFunc polynomial(x, y, 5);       // local degree-4 Neville polynomial
SplineInterpRealFunc spline(x, y);               // natural cubic spline
BarycentricRationalInterp rational(x, y, 3);     // pole-free rational interpolant

Real linearValue = linear(query);
Real polynomialValue = polynomial(query);
Real polynomialErrorEstimate = polynomial.getLastErrorEst();
Real splineValue = spline(query);
Real splineDerivative = spline.Derivative(query);
Real rationalValue = rational(query);
```

The demo prints the interpolated value from each method at the same query point. The polynomial
object exposes Neville's last signed error estimate, while the spline can differentiate and
integrate the constructed cubic pieces.

### 13.2 Interpolate a 2D function on a grid

**Problem:** you have a rectangular table of height measurements over x/y grid coordinates. Estimate
the height between grid nodes with both bilinear and bicubic spline interpolation.

**Setup**

```cpp
Vector<Real> gridX({ 0.0, 1.0, 2.0, 3.0, 4.0 });
Vector<Real> gridY({ 0.0, 0.8, 1.6, 2.4 });
Matrix<Real> z(5, 4, {
    10.0, 10.5, 11.0, 11.4,
    10.8, 11.4, 11.8, 12.1,
    11.3, 12.0, 12.2, 12.4,
    11.1, 11.7, 12.0, 12.2,
    10.6, 11.1, 11.5, 11.9
});
```

**Solve**

```cpp
BilinearInterp2D bilinear(gridX, gridY, z);
BicubicSplineInterp2D spline2D(gridX, gridY, z);

Real qx = 1.7, qy = 1.1;
Real bilinearValue = bilinear(qx, qy);
Real splineValue = spline2D(qx, qy);

Real zValue, dzdx, dzdy;
bilinear.interpWithDerivatives(qx, qy, zValue, dzdx, dzdy);
```

Bilinear interpolation is fast and only $C^0$ across cell boundaries. `BicubicSplineInterp2D`
builds splines along the grid directions and is smoother for tabulated smooth surfaces.

### 13.3 Interpolate a parametric curve

**Problem:** turn a short set of 2D waypoints into a callable parametric curve, using the same
normalized parameter interval $t\in[0,1]$ for a polyline and a smooth spline curve.

**Setup**

```cpp
Matrix<Real> points(5, 2, {
    0.0, 0.0,
    1.0, 0.2,
    1.8, 1.0,
    2.8, 0.8,
    3.5, 1.6
});
```

**Solve**

```cpp
LinInterpParametricCurve<2> polyline(points);
SplineInterpParametricCurve<2> smooth(points);

VectorN<Real, 2> polylinePoint = polyline(0.5);
VectorN<Real, 2> smoothPoint = smooth(0.5);
```

The linear curve preserves the waypoint polyline exactly. The spline curve gives a smoother path
through the same points by fitting each coordinate against chord-length parameter.

**Notes**

- Use low local polynomial order for equispaced real-function data; high-order global-looking
  interpolation is where Runge oscillations show up.
- Use splines for smooth tabulated data when local control and stable behavior matter more than
  a single high-degree polynomial.
- 2D interpolation expects a rectangular grid: `z(i,j)` stores the value at
  `(gridX[i], gridY[j])`.
- Runnable version: `Cookbook_Recipe13_InterpolatingFunctions()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 14: Curve fitting, including fitting to points

The curve-fitting API is in `<mml/algorithms/CurveFitting.h>`. This chapter covers ordinary
least squares, weighted diagonal least squares, general linear least squares with custom basis
functions, polynomial fitting, and ridge/Tikhonov regularization.

Important scope note: MML does **not** yet expose full generalized least squares with a full
observation covariance matrix for correlated errors. That missing capability is tracked separately
as Beads issue `MinimalMathLibrary-m7up`; the examples below use the APIs that exist today.

### 14.1 Fit a sensor calibration line

**Problem:** a thermistor circuit reports voltage, and a reference thermometer reports temperature.
Fit a calibration line and predict a temperature from a new voltage reading.

**Setup**

```cpp
#include <mml/algorithms/CurveFitting.h>

Vector<Real> voltage({ 0.12, 0.38, 0.64, 0.91, 1.17, 1.44, 1.70, 1.96, 2.23, 2.49 });
Vector<Real> temperatureC({ 3.1, 9.8, 16.4, 23.2, 29.9, 36.7, 43.1, 49.8, 56.6, 63.0 });
```

**Fit**

```cpp
auto fit = LinearLeastSquaresDetailed(voltage, temperatureC);
Real predicted = fit.a * 1.85 + fit.b;
```

`fit.a` is the slope and `fit.b` is the intercept in $T=aV+b$. The detailed result also reports
the residual norm, mean squared error, and $R^2$.

### 14.2 Use weighted least squares for unequal measurement uncertainty

**Problem:** a load-cell calibration has one questionable high-load reading. Keep the measurement,
but downweight it instead of letting it dominate the fit.

**Setup**

```cpp
Vector<Real> referenceKg({ 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0 });
Vector<Real> sensorMv({ 0.03, 1.02, 2.03, 3.05, 4.10, 5.02, 6.38, 6.95 });
Vector<Real> weights({ 4.0, 4.0, 4.0, 4.0, 3.0, 3.0, 0.5, 2.0 });
Vector<std::function<Real(Real)>> basis({
    [](Real) { return 1.0; },
    [](Real x) { return x; }
});
```

**Fit**

```cpp
auto fit = WeightedGeneralLinearLeastSquares(referenceKg, sensorMv, weights, basis);
Real predicted = fit.evaluate(9.0, basis);
```

This minimizes $\sum_i w_i(y_i-c_0-c_1x_i)^2$. Use weights proportional to confidence or inverse
variance when that information is known.

### 14.3 Fit a polynomial performance curve

**Problem:** pump manufacturers often publish head versus flow-rate data. Fit a quadratic curve
so the head can be estimated between measured operating points.

**Setup**

```cpp
Vector<Real> flowLps({ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 });
Vector<Real> headM({ 39.9, 39.2, 38.1, 36.5, 34.7, 32.3, 29.5, 26.2, 22.7 });
```

**Fit and evaluate**

```cpp
auto fit = PolynomialFit(flowLps, headM, 2);
Real headAt55 = EvaluatePolynomial(5.5, fit.coefficients);
```

The polynomial coefficients are stored as `[c0,c1,c2]`, so the model is
$h(q)=c_0+c_1q+c_2q^2$. Inspect the condition number when raising the degree.

### 14.4 Fit seasonal demand with Fourier basis functions

**Problem:** monthly energy demand has an annual cycle. Fit a compact seasonal model using sine
and cosine terms.

**Setup**

```cpp
Vector<Real> monthAngle(12);
for (int month = 0; month < 12; ++month)
    monthAngle[month] = 2.0 * Constants::PI * month / 12.0;

Vector<Real> energyMWh({ 42.1, 39.8, 36.5, 32.4, 29.6, 28.3,
                         30.1, 33.7, 37.9, 41.8, 45.6, 46.9 });
Vector<std::function<Real(Real)>> basis({
    [](Real) { return 1.0; },
    [](Real t) { return std::sin(t); },
    [](Real t) { return std::cos(t); },
    [](Real t) { return std::sin(2.0 * t); },
    [](Real t) { return std::cos(2.0 * t); }
});
```

**Fit**

```cpp
auto fit = GeneralLinearLeastSquares(monthAngle, energyMWh, basis);
Real seasonalAmplitude = std::sqrt(fit.coefficients[1] * fit.coefficients[1]
                                + fit.coefficients[2] * fit.coefficients[2]);
Real julyForecast = fit.evaluate(Constants::PI, basis);
```

The model is nonlinear-looking in time, but it is linear in the unknown coefficients, so ordinary
general linear least squares applies.

### 14.5 Fit an exponential response curve

**Problem:** a heated object cools toward ambient temperature. Fit a two-timescale response using
fixed exponential basis functions.

**Setup**

```cpp
Vector<Real> minutes({ 0.0, 2.0, 4.0, 7.0, 10.0, 14.0, 18.0, 24.0, 30.0, 40.0 });
Vector<Real> temperatureC({ 92.0, 78.6, 68.8, 57.1, 49.6, 42.9, 38.5, 34.4, 31.9, 29.6 });
Vector<std::function<Real(Real)>> basis({
    [](Real) { return 1.0; },
    [](Real t) { return std::exp(-t / 5.0); },
    [](Real t) { return std::exp(-t / 20.0); }
});
```

**Fit**

```cpp
auto fit = GeneralLinearLeastSquares(minutes, temperatureC, basis);
Real temperatureAt12 = fit.evaluate(12.0, basis);
```

Because the decay constants are chosen up front, this is still a linear least-squares problem in
the coefficients. Estimating the decay constants themselves would be nonlinear optimization.

### 14.6 Fit a smooth 2D path from points

**Problem:** GPS or robot-localization points are noisy. Fit a smooth path by parameterizing the
points by chord length, then fitting $x(t)$ and $y(t)$ separately.

**Setup**

```cpp
Vector<Real> x({ 0.0, 1.2, 2.5, 3.9, 5.1, 6.4, 7.2, 8.0 });
Vector<Real> y({ 0.0, 0.4, 1.1, 1.5, 1.2, 0.6, -0.1, -0.4 });
Vector<Real> t(x.size());
t[0] = -1.0;
Real totalLength = 0.0;
for (int i = 1; i < x.size(); ++i) {
    Real dx = x[i] - x[i - 1];
    Real dy = y[i] - y[i - 1];
    totalLength += std::sqrt(dx * dx + dy * dy);
    t[i] = totalLength;
}
for (int i = 1; i < t.size(); ++i)
    t[i] = -1.0 + 2.0 * t[i] / totalLength;
```

**Fit**

```cpp
auto basis = MakeLegendreFitBasis(3);
auto fitX = RidgeGeneralLinearLeastSquares(t, x, basis, Real{1e-4});
auto fitY = RidgeGeneralLinearLeastSquares(t, y, basis, Real{1e-4});

Real midX = fitX.evaluate(0.0, basis);
Real midY = fitY.evaluate(0.0, basis);
```

Legendre basis functions are natural on $[-1,1]$, so the chord-length parameter is scaled to that
interval. Ridge regularization damps unnecessary wiggle in the fitted path.

**Notes**

- Curve fitting approximates noisy data; interpolation is for passing exactly through data.
- Use `LinearLeastSquaresDetailed` for simple calibration lines and diagnostics.
- Use `GeneralLinearLeastSquares` when the model is linear in coefficients but uses custom basis
  functions.
- Use weighted least squares for independent observations with unequal confidence; true correlated
  covariance GLS is tracked separately as `MinimalMathLibrary-m7up`.
- Runnable version: `Cookbook_Recipe14_CurveFitting()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 15: Curves and surfaces

MML has two layers for curves and surfaces. Use the ready-made classes in `<mml/core/Curves.h>`
and `<mml/core/Surfaces.h>` when the shape is one of the standard examples. Use the generic
parametric interfaces, or derive from `Surfaces::ISurfaceCartesian`, when the geometry is your own.

This recipe focuses on geometry objects as mathematical objects: sample them, differentiate them,
measure length and curvature, build Frenet frames for curves, and inspect first and second
fundamental forms for surfaces.

### 15.1 Use built-in curves and surfaces

**Problem:** evaluate a few common curves and surfaces without writing the parametrization by hand.

**Setup**

```cpp
#include <mml/core/Curves.h>
#include <mml/core/Surfaces.h>

Curves::Circle2DCurve circle(2.0, Pnt2Cart(1.0, -1.0));
Curves::HelixCurve helix(1.5, 0.4);
Surfaces::Sphere sphere(3.0);
Surfaces::Torus torus(3.0, 0.75);
```

**Solve**

```cpp
VectorN<Real, 2> circlePoint = circle(Constants::PI / 3.0);
VectorN<Real, 3> helixPoint = helix(Constants::PI);
VectorN<Real, 3> spherePoint = sphere(Constants::PI / 2.0, Constants::PI / 4.0);
VectorN<Real, 3> torusPoint = torus(0.0, Constants::PI);
```

`Circle2DCurve` and `HelixCurve` are curves: they take one parameter. `Sphere`, `Torus`,
`CylinderSurface`, `Helicoid`, `EnneperSurface`, and the other surface classes take two parameters.
The parameter names are conventionally `u` and `w`; check each class's domain before sampling near
singularities such as sphere poles.

### 15.2 Define custom curves and surfaces

**Problem:** define your own trajectory and surface while still using MML's derivative and geometry
helpers.

**Setup**

```cpp
VectorN<Real, 3> CustomCurvePoint(Real t)
{
    return VectorN<Real, 3>{ t, std::sin(t), Real{0.25} * t * t };
}

class WavySheet : public Surfaces::ISurfaceCartesian
{
public:
    Real getMinU() const override { return -2.0; }
    Real getMaxU() const override { return 2.0; }
    Real getMinW() const override { return -2.0; }
    Real getMaxW() const override { return 2.0; }

    VectorN<Real, 3> operator()(Real u, Real w) const override
    {
        return VectorN<Real, 3>{ u, w, Real{0.20} * std::sin(u) * std::cos(w) };
    }
};

Curves::CurveCartesian3D trajectory(0.0, 4.0, CustomCurvePoint);
WavySheet sheet;
```

**Solve**

```cpp
VectorN<Real, 3> point = trajectory(2.0);
VectorN<Real, 3> tangent = trajectory.getTangent(2.0);
VectorN<Real, 3> sheetPoint = sheet(1.0, 0.5);
VectorN<Real, 3> sheetNormal = sheet.Normal(1.0, 0.5);
```

For quick curves, `Curves::CurveCartesian3D` is enough. For surfaces where you want normals,
curvatures, principal directions, and fundamental forms, derive from `Surfaces::ISurfaceCartesian`.

### 15.3 Calculate basic curve properties

**Problem:** measure curvature, torsion, and arc length from curve objects.

**Setup**

```cpp
#include <mml/core/Integration/PathIntegration.h>

Curves::Circle2DCurve circle(2.0);
Curves::HelixCurve helix(1.5, 0.4);
```

**Solve**

```cpp
Real circleCurvature = circle.getCurvature(Constants::PI / 4.0);
Real helixLength = PathIntegration::ParametricCurveLength<3>(helix, 0.0, 2.0 * Constants::PI);
Real helixCurvature = helix.getCurvature(1.0);
Real helixTorsion = helix.getTorsion(1.0);
```

For a circle of radius 2, the curvature is $1/R=0.5$. For a helix, curvature and torsion are
constant when radius and pitch are constant, while arc length is computed by numerical integration
over the requested parameter interval.

### 15.4 Compute a Frenet frame

**Problem:** compute the moving trihedron of a 3D curve: tangent $T$, principal normal $N$, and
binormal $B$.

**Setup**

```cpp
#include <mml/core/DifferentialGeometry/Frames.h>

Curves::HelixCurve helix(1.5, 0.4);
Real t = 1.0;
```

**Solve**

```cpp
Vector3Cartesian tangent, normal, binormal;
helix.getMovingTrihedron(t, tangent, normal, binormal);

auto frame = DifferentialGeometry::ComputeFrenetFrame<3>(helix, t);
Real speed = frame.speed;
Real curvature = frame.curvature;
```

`getMovingTrihedron` is the curve-class convenience API. `ComputeFrenetFrame` works from the
generic `IParametricCurve<N>` interface and also reports speed and curvature. Avoid points where
curve speed is zero; the Frenet frame is undefined at singular parameters.

### 15.5 Inspect surface normals and curvatures

**Problem:** classify and compare surfaces by normals, Gaussian curvature, mean curvature, and
principal curvatures.

**Setup**

```cpp
Surfaces::Sphere sphere(2.0);
Surfaces::CylinderSurface cylinder(1.5, 4.0);
Surfaces::Helicoid helicoid(0.4);

Real sphereU = Constants::PI / 3.0;
Real sphereW = Constants::PI / 5.0;
```

**Solve**

```cpp
Real k1, k2;
sphere.PrincipalCurvatures(sphereU, sphereW, k1, k2);

VectorN<Real, 3> normal = sphere.Normal(sphereU, sphereW);
Real gaussian = sphere.GaussianCurvature(sphereU, sphereW);
Real mean = sphere.MeanCurvature(sphereU, sphereW);
Real cylinderK = cylinder.GaussianCurvature(0.7, 2.0);
Real helicoidH = helicoid.MeanCurvature(1.0, 0.7);
```

Gaussian curvature is intrinsic. Mean curvature depends on the embedding and normal convention.
For the current MML surface convention, the sphere's Gaussian curvature is positive and its mean
curvature is negative with the normal orientation used by the implementation.

### 15.6 Calculate first and second fundamental forms

**Problem:** extract the metric and curvature coefficients of a parametric surface.

**Setup**

```cpp
Surfaces::Torus torus(3.0, 0.75);
Real u = 0.8;
Real w = 1.1;
```

**Solve**

```cpp
Real E, F, G;
Real L, M, N;
torus.GetFirstNormalFormCoefficients(u, w, E, F, G);
torus.GetSecondNormalFormCoefficients(u, w, L, M, N);

Real K = torus.GaussianCurvature(u, w);
Real H = torus.MeanCurvature(u, w);
```

The first fundamental form $I=(E,F,G)$ is the surface metric:
$ds^2=E\,du^2+2F\,du\,dw+G\,dw^2$. The second fundamental form $II=(L,M,N)$ measures how the
surface bends relative to its normal. MML computes both numerically from the surface
parametrization, so stay away from degenerate parameter points where the two tangent directions
become parallel.

**Notes**

- Built-in curve classes include circles, spirals, a helix, twisted cubic, toroidal spiral, and
  several plane curves useful for visualization and geometry tests.
- Built-in surface classes include planes, cylinders, spheres, tori, ellipsoids, hyperboloids,
  paraboloids, helicoids, Enneper surfaces, Dini surfaces, and non-orientable examples.
- Use fitting or interpolation recipes when the curve comes from data. Use this recipe when the
  curve or surface is already a mathematical parametrization.
- Runnable version: `Cookbook_Recipe15_CurvesAndSurfaces()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 16: Polynoms: polynomial models and approximations

`Polynom<CoefT, FieldT>` in `<mml/base/Polynom.h>` stores coefficients in ascending powers:
`{a0,a1,a2}` means $a_0+a_1x+a_2x^2$. This chapter treats polynomials as reusable numerical
objects: evaluate them, combine them, divide them, differentiate/integrate them, build them from
measurements, and convert stable Chebyshev approximations back to ordinary power-basis polynomials.

Polynomial root finding is covered earlier in [Recipe 2.3](#23-polynomial-roots); the examples
below focus on what else polynomial objects are good for.

### 16.1 Build and evaluate a polynomial model

**Problem:** represent a small calibration correction and evaluate it without losing track of
coefficient order.

**Setup**

```cpp
#include <mml/base/Polynom.h>

Polynom<Real> calibration({ 1.20, -0.35, 0.08, -0.004 });
const Real input = 8.0;
```

**Solve**

```cpp
Real corrected = calibration(input);
Real constant = calibration[0];
int degree = calibration.degree();
Real leading = calibration.leadingTerm();
// p(8) = 1.472
```

The initializer is not descending powers. The first element is the constant term, the second is
the coefficient of $x$, and the last element is the leading coefficient.

### 16.2 Compose polynomial arithmetic into a power model

**Problem:** combine simple load-dependent voltage, current, and loss models into one net-power
polynomial.

**Setup**

```cpp
Polynom<Real> voltage({ 12.0, -0.30 });       // 12 - 0.30x
Polynom<Real> current({ 1.50, 0.20, -0.01 }); // 1.50 + 0.20x - 0.01x^2
Polynom<Real> cableLoss({ 0.0, 0.0, 0.35 });  // 0.35x^2
const Real load = 4.0;
```

**Solve**

```cpp
Polynom<Real> netPower = voltage * current - cableLoss;
Real available = netPower(load);
// netPower has degree 3; netPower(4) = 17.512
```

Polynomial addition, subtraction, multiplication, and scalar operations make small model
composition explicit and easy to audit.

### 16.3 Divide polynomials to inspect quotient and residual

**Problem:** separate a measured cubic trend into a known first-order factor, a quotient trend,
and the residual left over.

**Setup**

```cpp
Polynom<Real> knownFactor({ -2.0, 1.0 });       // x - 2
Polynom<Real> processTrend({ 3.0, 0.5, -0.1 }); // 3 + 0.5x - 0.1x^2
Polynom<Real> measured = knownFactor * processTrend + Polynom<Real>::Constant(0.25);
```

**Solve**

```cpp
Polynom<Real> quotient, remainder;
Polynom<Real>::poldiv(measured, knownFactor, quotient, remainder);

Polynom<Real> reconstructed = knownFactor * quotient + remainder;
Real check = reconstructed(5.0);
// quotient degree = 2, remainder = 0.25, check = measured(5)
```

The division identity is $u(x)=v(x)q(x)+r(x)$, with `degree(r) < degree(v)` unless the remainder
is zero.

### 16.4 Differentiate and integrate polynomial models

**Problem:** use one height polynomial to compute height, velocity, acceleration, and accumulated
height over a time interval exactly in polynomial arithmetic.

**Setup**

```cpp
Polynom<Real> height({ 120.0, 35.0, -4.9, 0.15 }); // h(t)
const Real t = 3.0;
```

**Solve**

```cpp
Polynom<Real> velocity = height.derivative();
Polynom<Real> acceleration = velocity.derivative();
Polynom<Real> accumulatedHeight = height.integral();

Vector<Real> jet(4);
height.Derive(t, jet); // p, p', p'', p''' at t

Real h = jet[0];
Real v = velocity(t);
Real a = acceleration(t);
Real integral04 = accumulatedHeight(4.0) - accumulatedHeight(0.0);
```

`derivative()` and `integral()` return new polynomial objects. `Derive(x, jet)` is useful when
several derivatives are needed at the same point.

### 16.5 Construct an explicit polynomial from measurements

**Problem:** convert a small exact calibration table into one polynomial object, then reuse it for
evaluation and sensitivity.

**Setup**

```cpp
std::vector<Real> dose({ 0.0, 1.0, 2.0, 3.0, 4.0 });
std::vector<Real> response({ 2.0, 3.4, 4.0, 4.1, 4.4 });
```

**Solve**

```cpp
Polynom<Real> responseCurve = Polynom<Real>::FromValues(dose, response);
Polynom<Real> sensitivity = responseCurve.derivative();

Real predicted = responseCurve(2.5);
Real slope = sensitivity(2.5);
// degree = 4, p(2.5) = 4.078125, p'(2.5) = 0.0708333...
```

Use this when you specifically need the explicit polynomial object. For many noisy points,
Recipe 14's fitting tools or Recipe 13's splines are usually more appropriate.

### 16.6 Use Chebyshev approximation and convert to a polynomial

**Problem:** approximate a smooth function stably on $[-1,1]$, then convert the Chebyshev series
to an ordinary power-basis `Polynom<Real>` for APIs that expect one.

**Setup**

```cpp
#include <mml/base/ChebyshevApproximation.h>

auto target = [](Real x) { return std::exp(x); };
ChebyshevApproximation approximation(target, -1.0, 1.0, 12);
const Real x = 0.35;
```

**Solve**

```cpp
Polynom<Real> powerSeries = approximation.ToPolynomial();

Real exact = target(x);
Real chebError = std::abs(approximation(x) - exact);
Real powerError = std::abs(powerSeries(x) - exact);
Real maxError = approximation.MaxError(target, 101);
```

On `[-1,1]`, `ToPolynomial()` returns a polynomial in the same variable. For a general interval
`[a,b]`, the Chebyshev object evaluates in the original variable, while `ToPolynomial()` returns
the power series in the mapped variable $y=(2x-a-b)/(b-a)$.

**Notes**

- `Polynom` coefficient order is ascending powers: index equals exponent.
- Use polynomial division to check known factors or isolate residual behavior.
- `FromValues` interpolates exactly; it is not a noisy-data fitting routine.
- Chebyshev approximations are often a better construction route for smooth functions than
  directly manipulating high-degree power-basis coefficients.
- Runnable version: `Cookbook_Recipe16_Polynoms()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 17: Quaternions

MML uses Hamilton quaternions stored as $[w,x,y,z]$. Rotation quaternions are unit length,
angles are in radians, and `Rotate` performs an active right-handed rotation
$v'=q[0,v]q^{-1}$. The two quaternions $q$ and $-q$ represent the same orientation.

### 17.1 Rotate a vector around an arbitrary axis

**Problem:** rotate the vector $(0,0,1)$ by 90° around the diagonal axis $(1,1,0)$.

**Setup** — `FromAxisAngle` expects a unit axis:

```cpp
#include <mml/base/Quaternions.h>

Vec3Cart axis = Vec3Cart(1.0, 1.0, 0.0).Normalized();
Vec3Cart original(0.0, 0.0, 1.0);
```

**Solve**

```cpp
Quaternion rotation = Quaternion::FromAxisAngle(axis, Constants::PI / 2.0);
Vec3Cart rotated = rotation.Rotate(original);
// rotated = (sqrt(2)/2, -sqrt(2)/2, 0); its length remains 1
```

Quaternion rotation preserves vector length when the quaternion is unit length.

### 17.2 Point one direction toward another

**Problem:** construct the shortest rotation that turns an object's forward direction toward a
target direction. The cross product supplies the axis and the dot product supplies the angle.

```cpp
Vec3Cart from = fromDirection.Normalized();
Vec3Cart to = toDirection.Normalized();
Real dot = std::clamp(from[0] * to[0] + from[1] * to[1] + from[2] * to[2], -1.0, 1.0);

Quaternion aim;
if (dot > 1.0 - 1e-12) {
    aim = Quaternion::Identity();                         // already aligned
} else if (dot < -1.0 + 1e-12) {
    Vec3Cart helper = std::abs(from[0]) < 0.9 ? Vec3Cart(1,0,0) : Vec3Cart(0,1,0);
    aim = Quaternion::FromAxisAngle(VectorProduct(from, helper).Normalized(), Constants::PI);
} else {
    Vec3Cart axis = VectorProduct(from, to).Normalized();
    aim = Quaternion::FromAxisAngle(axis, std::acos(dot));
}

Vec3Cart aligned = aim.Rotate(from); // equals to, within numerical precision
```

Parallel vectors have no unique axis but need no rotation. Antiparallel vectors need 180° around
*any* perpendicular axis, so that case must not normalize the zero cross product.

### 17.3 Compose rotations in the correct order

**Problem:** apply a 90° yaw around Z, then a 90° roll around X. Quaternion multiplication is
noncommutative; for active rotations the rightmost quaternion acts first.

```cpp
Quaternion yawZ = Quaternion::FromAxisAngle(Vec3Cart(0,0,1), Constants::PI / 2.0);
Quaternion rollX = Quaternion::FromAxisAngle(Vec3Cart(1,0,0), Constants::PI / 2.0);
Vec3Cart vector(1,0,0);

Vec3Cart yawThenRoll = (rollX * yawZ).Rotate(vector); // approximately (0,0,1)
Vec3Cart rollThenYaw = (yawZ * rollX).Rotate(vector); // approximately (0,1,0)
```

In general, to apply `first` and then `second`, compose `second * first`.

### 17.4 Convert between Euler angles, quaternions, and matrices

**Problem:** accept yaw-pitch-roll input, use a quaternion internally, and export a rotation matrix
for an API that expects $R v$ with column vectors.

```cpp
const Real deg = Constants::PI / 180.0;
Quaternion orientation = Quaternion::FromEulerZYX(30 * deg, 20 * deg, 10 * deg);

Vec3Cart yawPitchRoll = orientation.ToEulerZYX();
MatrixNM<Real, 3, 3> matrix = orientation.ToRotationMatrix();
Quaternion reconstructed = Quaternion::FromRotationMatrix(matrix);
// yawPitchRoll / deg = (30, 20, 10)
```

`FromEulerZYX(yaw,pitch,roll)` means $R=R_zR_yR_x$. Euler extraction is not unique and becomes
singular at pitch ±90°; retain the quaternion as the authoritative orientation.

### 17.5 Undo rotations and calculate relative orientation

**Problem:** undo an orientation and find the rotation that carries one orientation into another.

```cpp
const Real deg = Constants::PI / 180.0;
Quaternion from = Quaternion::FromAxisAngle(Vec3Cart(0,0,1), 30 * deg);
Quaternion to = Quaternion::FromAxisAngle(Vec3Cart(0,0,1), 100 * deg);

Quaternion relative = to * from.Inverse(); // 70 degrees around Z

Vec3Cart vector(2,-1,3);
Vec3Cart restored = from.Inverse().Rotate(from.Rotate(vector)); // original vector
```

For unit quaternions `Inverse()` equals `Conjugate()`, but `Inverse()` also works for general
nonzero quaternions by dividing by the squared norm.

### 17.6 Interpolate orientations with SLERP

**Problem:** smoothly animate from no rotation to 120° around Z at constant angular speed.

```cpp
const Real deg = Constants::PI / 180.0;
Quaternion start = Quaternion::Identity();
Quaternion finish = Quaternion::FromAxisAngle(Vec3Cart(0,0,1), 120 * deg);

Quaternion quarter = Quaternion::Slerp(start, finish, 0.25); // 30 degrees
Quaternion halfway = Quaternion::Slerp(start, finish, 0.50); // 60 degrees

Quaternion lerpQuarter = Quaternion::Lerp(start, finish, 0.25).Normalized();
// normalized LERP gives about 27.80 degrees, not constant angular speed
```

`Slerp` returns unit quaternions and negates the second endpoint when needed to follow the shortest
orientation path. Raw `Lerp` is cheaper but must be normalized before rotation use.

**Notes**

- Normalize externally supplied quaternion components before calling `Rotate`.
- Quaternion components use `[w,x,y,z]`, not the `[x,y,z,w]` order common in some graphics APIs.
- MML uses Hamilton multiplication; libraries using the JPL convention have different signs.
- Compare orientations through their rotated vectors or accept both $q$ and $-q$; component-wise
  equality alone does not capture orientation equivalence.
- Runnable version: `Cookbook_Recipe17_Quaternions()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 18: Analyzing functions

`RealFunctionAnalyzer` in `<mml/algorithms/Analyzers/FunctionsAnalyzer.h>` combines sampling,
numerical differentiation, root finding, limits, integration, and statistics to diagnose real
functions. Its conclusions are numerical evidence over a chosen interval and tolerance, not
symbolic proofs.

### 18.1 Produce a point and interval health report

**Problem:** inspect the domain, continuity, derivative, monotonicity, and sampled range of
$f(x)=\log x$, while safely detecting overflow in another function.

```cpp
#include <mml/algorithms/Analyzers/FunctionsAnalyzer.h>

RealFunctionFromStdFunc logarithm([](Real x) { return std::log(x); });
RealFunctionAnalyzer analyzer(logarithm, "log(x)");

auto valid = analyzer.AnalyzePoint(1.0);
auto invalid = analyzer.AnalyzePoint(-1.0);
auto interval = analyzer.AnalyzeInterval(0.25, 4.0, 200);
Real minimum = analyzer.MinInNPoints(0.25, 4.0, 200);
Real maximum = analyzer.MaxInNPoints(0.25, 4.0, 200);

RealFunctionFromStdFunc exponential([](Real x) { return std::exp(x); });
SafeEvalResult overflow = SafeEvaluate(exponential, 1000.0);
```

At $x=1$, log is defined, continuous, and differentiable; at $x=-1$ it is undefined. On
$[0.25,4]$ it is monotonic with sampled range $[-1.38629,1.38629]$. `SafeEvaluate` marks
$e^{1000}$ as overflow instead of treating infinity as a valid function value.

`AnalyzeInterval` samples `numPoints` positions using $x_i=x_1+i(x_2-x_1)/N$ for
$i=0,\ldots,N-1$, so it does not include the right endpoint. `MinInNPoints` and `MaxInNPoints`
sample both endpoints.

### 18.2 Build a complete polynomial feature map

**Problem:** locate the roots, extrema, inflection point, and monotonic regions of
$f(x)=x^3-3x$.

```cpp
RealFunctionFromStdFunc polynomial([](Real x) { return x*x*x - 3*x; });
RealFunctionAnalyzer analyzer(polynomial, "x^3-3x");

auto roots = analyzer.GetRoots(-2.0, 2.0, 1e-7);
auto extrema = analyzer.GetLocalOptimumsClassified(-2.0, 2.0, 1e-4);
auto inflections = analyzer.GetInflectionPoints(-2.0, 2.0, 1e-4);

bool leftBranch = analyzer.isMonotonic(-2.0, -1.0, 100); // true
bool wholeRange = analyzer.isMonotonic(-2.0, 2.0, 100);  // false
```

The analyzer finds roots $-\sqrt3,0,\sqrt3$, extrema at $-1$ and $1$, and an inflection at $0$.
Classified extrema include position, value, and `LOCAL_MINIMUM`/`LOCAL_MAXIMUM` type.

Root discovery scans 1000 intervals for near-zero samples or sign changes. Even-multiplicity roots
do not change sign and can be missed unless the sampling grid lands sufficiently near them.

### 18.3 Detect and classify discontinuities

**Problem:** distinguish jump, removable, and infinite discontinuities numerically.

```cpp
RealFunctionFromStdFunc step([](Real x) { return x < 0 ? 0.0 : 1.0; });
RealFunctionFromStdFunc removable([](Real x) {
    if (x == 1.0) return std::numeric_limits<Real>::quiet_NaN();
    return (x*x - 1.0) / (x - 1.0);
});
RealFunctionFromStdFunc pole([](Real x) { return 1.0 / (x - 2.0); });

RealFunctionAnalyzer stepAnalyzer(step);
RealFunctionAnalyzer removableAnalyzer(removable);
RealFunctionAnalyzer poleAnalyzer(pole);

auto jump = stepAnalyzer.ClassifyDiscontinuity(0.0);
auto hole = removableAnalyzer.ClassifyDiscontinuity(1.0);
auto infinite = poleAnalyzer.ClassifyDiscontinuity(2.0, 1e-8);
auto discovered = stepAnalyzer.FindDiscontinuities(-1.0, 1.0, 100);
```

The step has limits 0 and 1 and jump size 1. The removable hole has matching limits near 2 but an
undefined point value. The pole is classified `INFINITE`. `ComputeLeftLimit` and
`ComputeRightLimit` are also available directly.

Detection and classification are sampling/tolerance heuristics. Difficult oscillatory behavior,
very narrow features, or poles between coarse sample points may need a denser scan and tuned `eps`.

### 18.4 Estimate oscillation period from zero crossings

**Problem:** estimate the period of a damped oscillation without peak detection.

```cpp
RealFunctionFromStdFunc damped([](Real t) {
    return std::exp(-0.05*t) * std::sin(5*t);
});
RealFunctionAnalyzer analyzer(damped, "exp(-0.05t)sin(5t)");

auto roots = analyzer.GetRoots(0.2, 5.5, 1e-7); // 8 roots
Real zeroSpacing = analyzer.calcRootsPeriod(0.2, 5.5, 1000);
Real fullPeriod = 2.0 * zeroSpacing;
```

Average consecutive-zero spacing is 0.6283185, giving full period 1.2566371, matching $2\pi/5$.
Despite its name, `calcRootsPeriod` returns average **root spacing**. For a sine-like waveform this
is half a period; other waveforms can have a different number of zero crossings per cycle.

### 18.5 Measure approximation error

**Problem:** quantify how well a coarse three-node linear interpolation represents the parabola
$f(x)=9-(x-3)^2$ on $[0,6]$.

```cpp
RealFunctionFromStdFunc exact([](Real x) { return 9.0 - (x-3.0)*(x-3.0); });
Vector<Real> nodes({0.0, 3.0, 6.0});
Vector<Real> values({exact(0.0), exact(3.0), exact(6.0)});
LinearInterpRealFunc approximation(nodes, values);
RealFunctionComparer comparer(approximation, exact);

Real avg = comparer.getAbsDiffAvg(0.0, 6.0, 1000);       // about 1.5
Real maximum = comparer.getAbsDiffMax(0.0, 6.0, 1000);   // 2.25
Real signedArea = comparer.getIntegratedDiff(
    0.0, 6.0, IntegrationMethod::ROMBERG);                // -9
Real absoluteArea = comparer.getIntegratedAbsDiff(
    0.0, 6.0, IntegrationMethod::ROMBERG);                // 9
Real squaredError = comparer.getIntegratedSqrDiff(
    0.0, 6.0, IntegrationMethod::ROMBERG);                // 16.2
```

Signed error reveals bias; absolute error measures total deviation; squared error penalizes larger
local misses. Relative metrics are also available, but skip sample points where the reference is zero.

**Notes**

- Sampling resolution and tolerance are part of every result; rerun at finer resolution before
  trusting a narrow or rapidly varying feature.
- Sampled min/max are not guaranteed global optimization results.
- `PrintPointAnalysis`, `PrintIntervalAnalysis`, and `PrintDetailedIntervalAnalysis` provide
  ready-made human-readable reports over the same structured APIs.
- The header description mentions asymptotes and periodicity broadly, but there is no implemented
  high-level asymptote or symmetry detector; period support is root-spacing based.
- Runnable version: `Cookbook_Recipe18_AnalyzingFunctions()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 19: Coordinate transformations

> **The "coordinates to curvature" thread, part 1 of 4** (19 → 20 → 21 → 22): *where* you
> are → *what* acts there → *how* to measure → *what "straight" means* in a curved world.

Positions, velocities, and gradients all "change coordinates" — but **not by the same rule**.
This recipe establishes the three rules; Recipes 20 and 21 use them.

### 19.1 Points: Cartesian ↔ spherical ↔ cylindrical

The predefined transformation objects convert positions both ways (`transf` /
`transfInverse`); MML's spherical convention is (r, θ, φ) with θ the polar angle:

```cpp
#include <mml/core/CoordTransf/CoordTransfSpherical.h>
#include <mml/core/CoordTransf/CoordTransfCylindrical.h>

Vector3Cartesian p{ 1.0, 2.0, 2.0 };                       // |p| = 3 exactly

Vector3Spherical   pSpher = CoordTransfCartToSpher.transf(p);  // r = 3, theta, phi
Vector3Cylindrical pCyl   = CoordTransfCartToCyl.transf(p);

Vector3Cartesian backS = CoordTransfSpherToCart.transf(pSpher);   // round-trips exactly
```

### 19.2 Vectors are not points: a velocity transforms *contravariantly*

A velocity attached to a point does not transform like the point — its components follow
the Jacobian of the transformation, **evaluated at that point**:

```cpp
Vector3Cartesian x_cart{ 1.0, 2.0, 2.0 };
Vector3Spherical x_spher = CoordTransfCartToSpher.transf(x_cart);

Vec3Cart v_cart{ 1.0, 1.0, 0.0 };                          // |v|^2 = 2

Vector3Spherical v_spher = CoordTransfCartToSpher.transfVecContravariant(v_cart, x_cart);
Vector3Cartesian v_back  = CoordTransfSpherToCart.transfVecContravariant(v_spher, x_spher);
// v_back = (1, 1, 0) - round-trip to machine precision
```

The spherical components (1, 0.298, −0.2) look nothing like (1, 1, 0) — yet they describe
the same arrow. *What gives their numbers meaning? The metric — Recipe 21.*

### 19.3 Gradients transform by the other rule: *covariantly*

A gradient is not a velocity — it transforms with the **inverse** Jacobian. Using the 1/r
potential computed natively in spherical coordinates:

```cpp
#include <mml/core/Fields/Fields.h>
#include <mml/core/Fields/ScalarFieldOperations.h>

Vector3Cartesian p_cart{ 1.0, 2.0, 2.0 };
Vector3Spherical p_spher = CoordTransfCartToSpher.transf(p_cart);

ScalarFunction<3> potSpher(Fields::InverseRadialPotentialFieldSpher);
Vector3Spherical grad_spher = ScalarFieldOperations::GradientSpher(potSpher, p_spher);
// = (-1/r^2, 0, 0) = (-1/9, 0, 0)

Vector3Cartesian grad_cart = CoordTransfSpherToCart.transfVecCovariant(grad_spher, p_cart);
Vector3Spherical grad_back = CoordTransfCartToSpher.transfVecCovariant(grad_cart, p_spher);
// grad_back = (-1/9, 0, 0) - round-trip exact
```

**Notes**

- `transfVecContravariant` = velocity/displacement rule; `transfVecCovariant` =
  gradient/one-form rule. Mixing them up silently produces wrong physics.
- *Why two rules, and how are they related? The metric connects them — Recipe 21.3.*
- Runnable version: `Cookbook_Recipe19_CoordinateTransformations()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 20: Fields

> **The "coordinates to curvature" thread, part 2 of 4** — the physics must not care which
> coordinates you use. One field, every coordinate system, one answer.

### 20.1 One field, two coordinate systems, one physical answer

The gravitational/electrostatic potential φ = 1/r, differentiated independently in Cartesian
and in spherical coordinates. The covariant transform from Recipe 19.3 must reconcile them:

```cpp
Vector3Cartesian p_cart{ 1.0, 1.0, 1.0 };
Vector3Spherical p_spher = CoordTransfCartToSpher.transf(p_cart);

ScalarFunction<3> potCart(Fields::InverseRadialPotentialFieldCart);
ScalarFunction<3> potSpher(Fields::InverseRadialPotentialFieldSpher);

Vector3Cartesian grad_cart  = ScalarFieldOperations::GradientCart<3>(potCart, p_cart);
Vector3Spherical grad_spher = ScalarFieldOperations::GradientSpher(potSpher, p_spher);

Vector3Cartesian grad_transf = CoordTransfSpherToCart.transfVecCovariant(grad_spher, p_cart);
// grad_cart = grad_transf = (-0.19245, -0.19245, -0.19245): same physics, verified
```

### 20.2 The identities every field obeys

curl(grad φ) = 0 and div(curl A) = 0 — numerically, for non-trivial fields:

```cpp
#include <mml/core/Fields/VectorFieldOperations.h>

Vector3Cartesian p{ 1.2, -0.7, 0.4 };

ScalarFunction<3> phi([](const VectorN<Real, 3>& x) {
    return x[0] * x[0] * x[1] + x[2] * x[2] * x[2];
});
VectorFunctionFromStdFunc<3> gradPhi(std::function<VectorN<Real, 3>(const VectorN<Real, 3>&)>(
    [&phi](const VectorN<Real, 3>& x) {
        return ScalarFieldOperations::GradientCart<3>(phi, x);
    }));
Vec3Cart curlOfGrad = VectorFieldOperations::CurlCart(gradPhi, p);    // ~ (0, 0, 7e-14)

VectorFunction<3> A([](const VectorN<Real, 3>& x) {
    return VectorN<Real, 3>{ x[0] * x[0] * x[1], x[1] * x[1] * x[2], x[2] * x[2] * x[0] };
});
VectorFunctionFromStdFunc<3> curlA(std::function<VectorN<Real, 3>(const VectorN<Real, 3>&)>(
    [&A](const VectorN<Real, 3>& x) {
        return VectorN<Real, 3>(VectorFieldOperations::CurlCart(A, x));
    }));
Real divOfCurl = VectorFieldOperations::DivCart<3>(curlA, p);         // 0
```

### 20.3 Empty-space gravity: ∇²(1/r) = 0

The 1/r potential is *harmonic* — it satisfies the Laplace equation away from the source,
and the result must not depend on the coordinate system used to compute it:

```cpp
Real lapCart  = ScalarFieldOperations::LaplacianCart<3>(potCart, p_cart);    // ~ -8e-11
Real lapSpher = ScalarFieldOperations::LaplacianSpher(potSpher, p_spher);    // ~ -1e-11
```

**Notes**

- Gradient/divergence/curl/Laplacian exist in Cartesian, spherical, and cylindrical variants
  (`...Cart` / `...Spher` / `...Cyl`); each applies the correct scale factors internally.
- Predefined physics fields live in `Fields::` (inverse-radial potential and force in all
  three coordinate systems).
- *The scale factors in those formulas are exactly the metric components — Recipe 21.*
- Runnable version: `Cookbook_Recipe20_Fields()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 21: Tensors

> **The "coordinates to curvature" thread, part 3 of 4** — the machinery that gives
> components meaning, resolves Recipe 19's two transformation rules, and ends with a twist.

### 21.1 Tensor2 basics: index variance, contraction, evaluation

```cpp
#include <mml/base/Tensor/Tensor2.h>

Tensor2<3> T(1, 1, { 2, 0, 1,        // mixed tensor: 1 covariant, 1 contravariant index
                     0, 3, 0,
                     0, 0, 4 });

T.NumCovar();                        // 1
T.Contract();                        // trace = 9 (contraction needs one index of each kind)

VectorN<Real, 3> u{ 1, 0, 0 }, w{ 0, 0, 1 };
Real val = T(u, w);                  // tensor as a multilinear machine: = 1
```

### 21.2 THE tensor: the metric — what makes components meaningful

Recipe 19.2 left a puzzle: the velocity's spherical components (1, 0.298, −0.2) look nothing
like (1, 1, 0). The **metric tensor** is what turns components into geometry — |v|² = gᵢⱼvⁱvʲ
gives the same physical length in every coordinate system:

```cpp
#include <mml/core/MetricTensor.h>

MetricTensorSpherical g;
Tensor2<3> g_ij = g(x_spher);        // diag(1, r^2, r^2 sin^2 theta) = diag(1, 9, 5)

Real lenSq = 0.0;
for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
        lenSq += g_ij(i, j) * v_spher[i] * v_spher[j];
// lenSq = 2.000000 - exactly |v_cart|^2
```

### 21.3 Raising an index: gradient (covariant) → force (contravariant)

Recipe 19's two transformation rules are connected by the metric: g^ij raises the covariant
gradient into the contravariant force vector,

```cpp
auto g_inv = g.GetContravariantMetric(p_spher);

VectorN<Real, 3> force;              // F^i = g^ij grad_j
for (int i = 0; i < 3; i++) {
    force[i] = 0.0;
    for (int j = 0; j < 3; j++)
        force[i] += g_inv(i, j) * grad_spher[j];
}

// contravariant components use the VELOCITY rule from 19.2; carried to Cartesian
// (where the metric is the identity) they coincide with the Cartesian gradient:
Vector3Cartesian force_cart = CoordTransfSpherToCart.transfVecContravariant(
    Vector3Spherical(force), p_spher);
// force_cart = grad_cart - verified to 1e-7
```

### 21.4 The twist: curvy coordinates, FLAT space

Spherical coordinates *feel* curved — their Christoffel symbols are nonzero. But curvature
is a property of **space**, not of coordinates, and the Riemann tensor knows the difference:

```cpp
MetricTensorSpherical g;
Vector3Spherical pos{ 2.0, Constants::PI / 3, Constants::PI / 4 };

g.GetChristoffelSymbolSecondKind(0, 1, 1, pos);   // Gamma^r_theta,theta = -r = -2
g.GetChristoffelSymbolSecondKind(1, 0, 1, pos);   // Gamma^theta_r,theta = 1/r = 0.5

Tensor4<3> riemann = g.GetRiemannCurvatureTensor(pos);
// max |R^i_jkl| ~ 1e-10: the space is FLAT - all that "curvature" was the coordinates
```

**Notes**

- `MetricTensorField<N>` also provides `GetRicciTensor`, `GetRicciScalar`,
  `GetEinsteinTensor`, covariant derivatives, and `MetricTensorFromCoordTransf` builds the
  metric automatically from any coordinate transformation.
- *If Riemann ≠ 0, no coordinates can flatten it — that's real curvature. Recipe 22.*
- Runnable version: `Cookbook_Recipe21_Tensors()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 22: Differential geometry and geodesics

> **The "coordinates to curvature" thread, finale** — Recipe 21 showed that spherical
> coordinates only *pretended* to be curved. Now meet a space that is **truly** curved — the
> surface of a sphere — and find out what "straight" and "parallel" mean when you live on it.

### 22.1 The sphere's induced metric and the Theorema Egregium

A 2D creature on a sphere of radius R can measure its world's metric — induced from the
embedding — and from the metric *alone* compute the Gaussian curvature. It never needs to
see the third dimension (Gauss's **Theorema Egregium**):

```cpp
#include <mml/algorithms/Geodesic.h>     // brings InducedMetric2D + geodesic integrators

Surfaces::Sphere sphere(2.0);            // R = 2, parametrized by (theta, phi)
DifferentialGeometry::InducedMetric2D g(sphere);

VectorN<Real, 2> pos{ Constants::PI / 3, 0.5 };    // theta = 60 deg
Tensor2<2> g_ij = g(pos);                // diag(R^2, R^2 sin^2 theta) = diag(4, 3)

Real K = DifferentialGeometry::GaussianCurvatureIntrinsic(g, pos);
// K = 0.24989  (exact 1/R^2 = 0.25) - real curvature, no coordinates can remove it
```

### 22.2 Geodesics: the great circles

A geodesic is the straightest possible path — d²qⁱ/dλ² = −Γⁱⱼₖ q̇ʲq̇ᵏ, integrated by
`IntegrateGeodesicFixedStep`. Launch from the equator of a unit sphere heading north-east at
unit speed; the path is a great circle, and after arc length 2π it must come home:

```cpp
Surfaces::Sphere sphere(1.0);
DifferentialGeometry::InducedMetric2D metric(sphere);

VectorN<Real, 2> pos0{ Constants::PI / 2, 0.0 };               // on the equator
VectorN<Real, 2> vel0{ -1.0 / std::sqrt(2.0), 1.0 / std::sqrt(2.0) };  // NE, |v|_g = 1

auto sol = IntegrateGeodesicFixedStep<2>(metric, pos0, vel0, 0.0, 2 * Constants::PI, 2000);

auto end = sol.getXValuesAtEnd();        // (theta, phi, vtheta, vphi)
// theta = 1.5707960 (pi/2), phi = 6.2831854 (2 pi): back at the start point
// |v|_g^2 = 1.0000001: geodesics preserve speed - this is why planes fly great circles
```

### 22.3 Finale — holonomy: parallel transport remembers curvature

Carry a vector around the 60° latitude circle, keeping it as parallel as the surface allows
(`IntegrateParallelTransportFixedStep`). On a flat sheet it would come home unchanged. On
the sphere it comes home **rotated** — by exactly the curvature enclosed (Gauss-Bonnet:
α = K × cap area = 2π(1 − cos 60°) = **π**):

```cpp
const Real theta0 = Constants::PI / 3;             // the 60-deg colatitude circle
ParametricCurveFromStdFunc<2> latitude(std::function<VectorN<Real, 2>(Real)>(
    [theta0](Real t) { return VectorN<Real, 2>{ theta0, t }; }));

VectorN<Real, 2> V0{ 1.0, 0.0 };                   // unit vector pointing south

auto sol = IntegrateParallelTransportFixedStep<2>(metric, latitude, V0, 0.0, 2 * Constants::PI, 2000);
auto Vend = sol.getXValuesAtEnd();
// Vend = (-1.0000000, -0.0000000): the vector came home ANTIPARALLEL

// measure the rotation angle with the metric inner product:
// cos(angle) = g(V0, Vend) / (|V0| |Vend|)  ->  angle = 3.14159262  (= pi, Gauss-Bonnet)
```

The vector was never rotated locally — every step kept it maximally parallel. The rotation
is the curvature of the enclosed cap, made visible. *This* is the difference between
Recipe 21.4's flat space in curvy coordinates and a genuinely curved world.

**Notes**

- `InducedMetric2D` works for any `Surfaces::ISurfaceCartesian` (torus, Möbius strip, ...);
  `IntegrateSurfaceGeodesicFixedStep` is the one-call convenience for surfaces.
- Avoid coordinate singularities of the parametrization (the sphere's poles): Christoffel
  symbols blow up there even though the surface itself is perfectly regular.
- The same `GeodesicEquationSystem` machinery works for any `MetricTensorField<N>` — including
  Lorentzian metrics (`MetricTensorMinkowski`), where geodesics are free-fall world lines.
- Runnable version: `Cookbook_Recipe22_DifferentialGeometry()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 23: Fourier algorithms

Core transforms live in `<mml/algorithms/Fourier/Fourier.h>`. Real FFT, spectrum helpers,
windows, and convolution have focused headers in the same directory. The forward convention is

$$X_k=\sum_{n=0}^{N-1}x_n e^{-2\pi i kn/N}.$$

Legacy forward transforms are unscaled; legacy inverse transforms divide by $N$.

### 23.1 FFT round trip and normalization

**Problem:** move a complex signal to frequency space and back without changing its values.

```cpp
#include <mml/algorithms/Fourier/Fourier.h>

Vector<Complex> signal(16);
for (int i = 0; i < signal.size(); ++i)
  signal[i] = Complex(std::sin(0.3 * i), 0.25 * std::cos(0.7 * i));

auto spectrum = Fourier::FFT::Forward(signal);
auto recovered = Fourier::FFT::Inverse(spectrum); // legacy inverse applies 1/N

auto unitarySpectrum = Fourier::FFT::Forward(
  signal, Fourier::TransformNormalization::Orthonormal);
auto unitaryRecovered = Fourier::FFT::Inverse(
  unitarySpectrum, Fourier::TransformNormalization::Orthonormal);
```

Both round trips reproduce the input to about $6\times10^{-16}$. `FFT::Transform(data,-1)` is
the raw in-place inverse and must be divided by $N$ manually; the `Inverse` wrapper handles scaling.

### 23.2 Detect frequencies and amplitudes

**Problem:** recover the two tones in a real signal sampled at 256 Hz.

```cpp
#include <mml/algorithms/Fourier/FourierRealFFT.h>
#include <mml/algorithms/Fourier/FourierSpectrum.h>

const int N = 256;
const Real sampleRate = 256.0;
Vector<Real> signal(N);
for (int n = 0; n < N; ++n) {
  Real t = n / sampleRate;
  signal[n] = std::sin(2 * Constants::PI * 20 * t)
        + 0.4 * std::sin(2 * Constants::PI * 60 * t);
}

auto spectrum = Fourier::RealFFT::Forward(signal);
auto frequencies = Fourier::FrequencyAxis(
  N, sampleRate, Fourier::FrequencyAxisLayout::OneSided);
Real amplitudeAtBin = 2.0 * std::abs(spectrum[20]) / N;
```

The dominant bins are exactly 20 Hz at amplitude 1 and 60 Hz at amplitude 0.4. `RealFFT` stores
only $N/2+1$ bins: DC, positive frequencies, and Nyquist. DC and Nyquist are not doubled when
converting magnitudes to one-sided amplitudes.

### 23.3 Control spectral leakage with windows

**Problem:** estimate a 20.5-bin tone whose frequency does not fit an integer FFT bin. A rectangular
record has a discontinuity at its periodic boundary, spreading energy across the spectrum.

```cpp
#include <mml/algorithms/Fourier/FourierWindowing.h>

auto rectangular = Fourier::Windows::Rectangular(N);
auto hann = Fourier::Windows::Hann(N);
auto windowed = Fourier::Windows::ApplyWindow(signal, hann);
auto spectrum = Fourier::RealFFT::Forward(windowed);
auto metrics = Fourier::Windows::Metrics(hann);

int peakBin = 1;
for (int bin = 2; bin < spectrum.size(); ++bin)
  if (std::abs(spectrum[bin]) > std::abs(spectrum[peakBin])) peakBin = bin;

Real correctedPeak = 2.0 * std::abs(spectrum[peakBin])
           / (N * metrics.coherent_gain);
```

In the runnable comparison, power outside ±2 bins falls from 8.09% (rectangular) to 0.052%
(Hann). The tradeoff is a wider main lobe. Coherent gain corrects peak amplitude; ENBW measures
noise bandwidth. Flat-top gives high amplitude accuracy but has ENBW about 3.78 bins.

### 23.4 Remove interference in the frequency domain

**Problem:** recover a 10 Hz signal contaminated by 40 Hz and 70 Hz interference.

```cpp
auto spectrum = Fourier::RealFFT::Forward(noisy);
auto frequencies = Fourier::FrequencyAxis(
  N, sampleRate, Fourier::FrequencyAxisLayout::OneSided);

for (int bin = 0; bin < spectrum.size(); ++bin)
  if (frequencies[bin] > 20.0)
    spectrum[bin] = Complex{};

Vector<Real> filtered = Fourier::RealFFT::Inverse(spectrum);
```

For this exact-bin demonstration, low-pass filtering reduces RMSE from 0.515 to about
$2\times10^{-15}$. Real measurements need tapered transitions: an abrupt spectral cutoff corresponds
to a long oscillatory impulse response and can ring around transients.

### 23.5 Smooth signals with convolution

**Problem:** apply a three-tap smoothing kernel and choose the desired output extent.

```cpp
#include <mml/algorithms/Fourier/FourierConvolution.h>

Vector<Real> signal({1,2,3,4,5,4,3,2});
Vector<Real> kernel({0.25,0.5,0.25});

auto full = Fourier::Convolve(signal, kernel, Fourier::ConvolutionMode::Full);
auto same = Fourier::Convolve(signal, kernel, Fourier::ConvolutionMode::Same);
auto valid = Fourier::Convolve(signal, kernel, Fourier::ConvolutionMode::Valid);

auto sameFFT = Fourier::Convolve(signal, kernel,
  Fourier::ConvolutionMode::Same, Fourier::ConvolutionMethod::FFT);
```

Output sizes are 10, 8, and 6. Direct and FFT methods agree to about $4\times10^{-16}$.
`Auto` chooses FFT only when the estimated direct work is large enough to justify transform overhead.

### 23.6 Estimate time delay from phase

**Problem:** estimate a circular sample delay from the phase slope between two signals.

```cpp
const int N = 128, actualDelay = 7;
Vector<Complex> original(N, Complex{}), delayed(N, Complex{});
original[0] = 1.0;
delayed[actualDelay] = 1.0;

auto X = Fourier::FFT::Forward(original);
auto Y = Fourier::FFT::Forward(delayed);
Vector<Complex> ratio(N / 2 + 1);
for (int k = 0; k < ratio.size(); ++k)
  ratio[k] = Y[k] * std::conj(X[k]);

auto phase = Fourier::UnwrapPhase(Fourier::Phase(ratio));
// Fit phase[k] = slope*k; delay = -slope*N/(2*pi)
```

The fitted result is exactly 7 samples. Reliable estimates require spectral bins with meaningful
magnitude; phase is unstable where either signal has nearly zero energy.

### 23.7 Transform arbitrary sample counts with Bluestein

**Problem:** transform a prime-length signal without padding it to a power of two.

```cpp
const int N = 127;
Vector<Complex> signal(N);
for (int n = 0; n < N; ++n)
  signal[n] = Complex(std::sin(0.17*n) + 0.3*std::cos(0.41*n), 0.1*std::sin(0.07*n));

auto fast = Fourier::FFT::Forward(signal); // automatically dispatches to Bluestein
auto reference = Fourier::DFT::Forward(signal);
auto recovered = Fourier::FFT::Inverse(fast);
```

For $N=127$, FFT and the $O(N^2)$ reference DFT differ by about $1.3\times10^{-12}$; round-trip
error is about $3.6\times10^{-14}$. Complex FFT accepts any nonempty size. `RealFFT` currently
requires a power-of-two real input.

### 23.8 Verify energy with Parseval's theorem

**Problem:** prove that the transform accounts for all signal energy. For a one-sided real spectrum,
interior bins represent both positive and negative frequencies.

```cpp
Real timeEnergy = 0.0;
for (Real value : signal) timeEnergy += value * value;

auto spectrum = Fourier::RealFFT::Forward(signal);
Real frequencyEnergy = std::norm(spectrum[0]) + std::norm(spectrum[N/2]);
for (int k = 1; k < N/2; ++k)
  frequencyEnergy += 2.0 * std::norm(spectrum[k]);
frequencyEnergy /= N;
```

The example obtains 80 in both domains, differing by roughly $1.3\times10^{-13}$. Do not double
DC or Nyquist. `Fourier::Power` returns squared magnitudes; it is not a sample-rate/window-normalized
physical power spectral density.

### 23.9 Compress smooth data with the DCT

**Problem:** approximate 64 smooth samples using only eight cosine coefficients.

```cpp
Vector<Real> signal(64);
for (int n = 0; n < signal.size(); ++n) {
  Real x = (n + 0.5) / signal.size();
  signal[n] = std::exp(-2*x) + 0.2*std::cos(2*Constants::PI*x);
}

auto coefficients = Fourier::DCT::ForwardII(signal);
for (int k = 8; k < coefficients.size(); ++k) coefficients[k] = 0.0;
auto reconstructed = Fourier::DCT::InverseII(coefficients);
```

Keeping 8 of 64 coefficients gives RMSE 0.00771. DCT-II/DCT-III compact smooth, even-extended
data well and are common in compression and Chebyshev spectral approximation.

### 23.10 Analyze zero-boundary modes with the DST

**Problem:** identify vibration modes of samples whose implicit endpoints are fixed at zero.

```cpp
const int N = 32;
Vector<Real> signal(N);
for (int n = 0; n < N; ++n) {
  signal[n] = std::sin(Constants::PI * 3 * (n+1) / (N+1))
        + 0.3 * std::sin(Constants::PI * 7 * (n+1) / (N+1));
}

auto coefficients = Fourier::DCT::ForwardDST(signal);
auto recovered = Fourier::DCT::InverseDST(coefficients);
```

The dominant modes are 3 and 7; reconstruction RMSE is about $1.8\times10^{-15}$. MML's DST-I
uses an orthonormal sine basis and is its own inverse, matching homogeneous Dirichlet boundaries.

**Notes**

- `FrequencyAxis` supports one-sided, two-sided, and shifted layouts; use `FFTShift` on spectrum
  values when requesting a shifted two-sided axis.
- Zero padding samples a spectrum more densely but does not create new physical frequency resolution.
- Current core does not expose STFT/spectrogram, Welch PSD, 2D FFT, or correlation facades; those
  package/archive APIs are intentionally absent from these recipes.
- Runnable version: `Cookbook_Recipe23_FourierAlgorithms()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 24: Graph algorithms

`Graph<V,E>` in `<mml/base/Graph.h>` stores optional vertex payloads and edge weights in an
adjacency list. The aggregate `<mml/algorithms/GraphAlg.h>` provides traversal, paths, structure,
spanning-tree, flow, cut, and matching algorithms. Vertices are addressed by zero-based indices.

### 24.1 Build, traverse, and inspect a graph

**Problem:** discover which machines are reachable in a network containing three disconnected
regions.

```cpp
#include <mml/algorithms/GraphAlg.h>

Graph<std::string> network(Graph<std::string>::Type::Undirected);
for (const char* name : { "A", "B", "C", "D", "E", "F" })
    network.addVertex(name);
network.addEdge(0, 1);
network.addEdge(1, 2);
network.addEdge(0, 2);
network.addEdge(3, 4);

auto bfs = BFS(network, 0);                   // visit order A, B, C
auto dfs = DFS(network, 0);                   // depth-first traversal tree
auto components = ConnectedComponents(network); // 3 components
```

BFS distances are edge counts and its parent array reconstructs minimum-hop paths. DFS is useful
for structural exploration. Both visit only the component reachable from the chosen start.

### 24.2 Find the best route

**Problem:** compare “fewest road segments” with “shortest total distance” between a depot and
harbor. An A* heuristic uses straight-line coordinates to guide the weighted search.

```cpp
Graph<std::string> roads(Graph<std::string>::Type::Undirected);
for (const char* name : { "Depot", "Museum", "Park", "Station", "Harbor" })
    roads.addVertex(name);
roads.addEdge(0, 1, 2.0); roads.addEdge(0, 2, 4.0);
roads.addEdge(1, 2, 1.0); roads.addEdge(1, 3, 5.0);
roads.addEdge(2, 3, 1.5); roads.addEdge(3, 4, 2.0);
roads.addEdge(2, 4, 6.0);

auto hops = ShortestPathUnweighted(roads, 0, 4); // Depot -> Park -> Harbor
auto route = DijkstraPath(roads, 0, 4);           // weight 6.5

std::vector<Vec3Cart> position{{0,0,0}, {2,0,0}, {2,1,0}, {3,2,0}, {5,2,0}};
auto heuristic = [&](size_t v, size_t target) {
    return (position[target] - position[v]).NormL2();
};
auto guided = AStarPath(roads, 0, 4, heuristic);
```

Dijkstra and A* require nonnegative weights. A* remains optimal only when its heuristic never
overestimates the remaining cost. Use Bellman–Ford for negative edges and negative-cycle detection.

### 24.3 Schedule dependent tasks and find the critical path

**Problem:** order tasks subject to prerequisites and identify the dependency chain that determines
the earliest completion time.

```cpp
Graph<std::string> project(Graph<std::string>::Type::Directed);
for (const char* task : { "Start", "Design", "Procure", "Build", "Test", "Deploy" })
    project.addVertex(task);
project.addEdge(0, 1, 3.0); project.addEdge(0, 2, 2.0);
project.addEdge(1, 3, 5.0); project.addEdge(2, 3, 4.0);
project.addEdge(3, 4, 2.0); project.addEdge(4, 5, 1.0);

auto order = TopologicalSort(project);
auto critical = LongestPathDAG(project, 0);
auto criticalPath = critical.pathTo(5); // Start, Design, Build, Test, Deploy
// critical.distance[5] = 11
```

Topological order is generally not unique. `TopologicalSort` and `LongestPathDAG` report cycle
failure instead of returning a misleading schedule; `FindCycle` extracts an example cycle.

### 24.4 Design the cheapest connected network

**Problem:** connect five sites using the minimum total cable length.

```cpp
Graph<std::string> sites(Graph<std::string>::Type::Undirected);
for (const char* site : { "A", "B", "C", "D", "E" }) sites.addVertex(site);
sites.addEdge(0,1,4); sites.addEdge(0,2,2); sites.addEdge(1,2,1);
sites.addEdge(1,3,5); sites.addEdge(2,3,8); sites.addEdge(2,4,10);
sites.addEdge(3,4,2);

auto kruskal = Kruskal(sites);
auto prim = Prim(sites);
// Both return 4 edges with total weight 10
```

An MST minimizes total network cost, not the route between every pair. A disconnected input returns
a minimum spanning forest with `isComplete == false` and `DisconnectedGraph` status.

### 24.5 Find single points of failure

**Problem:** identify routers and links whose removal disconnects a network.

```cpp
Graph<std::string> routers(Graph<std::string>::Type::Undirected);
for (const char* name : { "A", "B", "C", "D", "E", "F" }) routers.addVertex(name);
routers.addEdge(0,1); routers.addEdge(1,2); routers.addEdge(2,0); // triangle
routers.addEdge(1,3);                                           // bridge
routers.addEdge(3,4); routers.addEdge(4,5); routers.addEdge(5,3); // triangle

auto resilience = UndirectedConnectivity(routers);
// articulationPoints = {1,3}; bridges = {{1,3}}; 3 biconnected regions
```

Articulation points are critical vertices; bridges are critical links. Inside a biconnected
component no single vertex failure disconnects the remaining vertices.

### 24.6 Analyze cycles and strongly connected subsystems

**Problem:** find circular dependencies and groups of services that can all reach one another.

```cpp
Graph<std::string> services(Graph<std::string>::Type::Directed);
for (const char* name : { "API", "Auth", "Users", "Billing", "Email", "Audit" })
    services.addVertex(name);
services.addEdge(0,1); services.addEdge(1,2); services.addEdge(2,0);
services.addEdge(2,3); services.addEdge(3,4); services.addEdge(4,3);
services.addEdge(4,5);

auto cycle = FindCycle(services);                    // API -> Auth -> Users -> API
auto scc = StronglyConnectedComponents(services);    // 3 SCCs
auto condensation = CondensationDAG(services);       // acyclic component graph
```

Collapsing each strongly connected component produces a DAG, which is often the useful high-level
view of a dependency, call, or state-transition graph.

### 24.7 Calculate maximum throughput and the bottleneck cut

**Problem:** find maximum pipe throughput and the saturated boundary that limits it.

```cpp
Graph<std::string> pipes(Graph<std::string>::Type::Directed);
for (const char* name : { "Source", "A", "B", "C", "D", "Sink" })
    pipes.addVertex(name);
pipes.addEdge(0,1,16); pipes.addEdge(0,2,13); pipes.addEdge(1,2,10);
pipes.addEdge(2,1,4);  pipes.addEdge(1,3,12); pipes.addEdge(2,4,14);
pipes.addEdge(3,2,9);  pipes.addEdge(4,3,7);  pipes.addEdge(3,5,20);
pipes.addEdge(4,5,4);

auto flow = Dinic(pipes, 0, 5);
// flow.maxFlow = 23; cutEdges is a minimum-capacity source/sink cut
```

Flow algorithms interpret weights as nonnegative capacities and require a directed graph. Dinic is
the faster general choice; `EdmondsKarp` is a simpler alternative with the same result structure.

### 24.8 Assign jobs to workers

**Problem:** assign as many workers as possible to compatible jobs, with each worker and job used
at most once.

```cpp
Graph<std::string> assignments(Graph<std::string>::Type::Undirected);
for (const char* name : { "Ana", "Bo", "Cy", "Weld", "Inspect", "Pack" })
    assignments.addVertex(name);
assignments.addEdge(0,3); assignments.addEdge(0,4);
assignments.addEdge(1,3); assignments.addEdge(1,5);
assignments.addEdge(2,4); assignments.addEdge(2,5);

auto matching = HopcroftKarp(assignments, std::vector<size_t>{0,1,2});
// cardinality = 3: every worker receives a compatible job
```

The explicit partition removes ambiguity about which vertices represent workers. Every edge must
cross the bipartition; invalid same-side edges return `GraphTypeMismatch`.

**Notes**

- Undirected edges appear in both adjacency lists but `numEdges()` counts each edge once.
- Duplicate edges and disallowed self-loops throw; self-loops can be enabled at construction.
- Check `succeeded()`, `graphStatus`, and algorithm-specific flags such as `found`, `isDAG`, or
  `isComplete` before consuming a result.
- `toAdjacencyMatrix()`, `toLaplacianMatrix()`, and `toNormalizedLaplacian()` connect graph
  structure to MML's matrix and eigensolver tools for later spectral analysis.
- Runnable version: `Cookbook_Recipe24_GraphAlgorithms()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 25: DAE solvers

MML's DAE solvers handle semi-explicit index-1 systems in the form

$$
\frac{d\mathbf{x}}{dt}=\mathbf{f}(t,\mathbf{x},\mathbf{y}),\qquad
\mathbf{0}=\mathbf{g}(t,\mathbf{x},\mathbf{y})
$$

where `x` contains differential variables and `y` contains algebraic variables. The index-1
condition is practical and important: $\partial g/\partial y$ must be nonsingular so the
constraints can determine the algebraic variables locally.

### 25.1 Define and solve an index-1 DAE

**Problem:** solve a one-variable DAE with the constraint $x+y=1$:

$$
x'=-x+y,\qquad 0=x+y-1
$$

```cpp
#include <mml/algorithms/DAESolvers.h>

class LinearDAE : public IODESystemDAEWithJacobian {
public:
  int getDiffDim() const override { return 1; }
  int getAlgDim() const override { return 1; }

  void diffEqs(Real, const Vector<Real>& x, const Vector<Real>& y,
               Vector<Real>& dxdt) const override {
    dxdt[0] = -x[0] + y[0];
  }

  void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y,
                      Vector<Real>& g) const override {
    g[0] = x[0] + y[0] - 1.0;
  }

  void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0,0) = -1.0; }
  void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0,0) =  1.0; }
  void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0,0) =  1.0; }
  void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0,0) =  1.0; }
};
```

```cpp
LinearDAE system;
Vector<Real> x0{1.0};
Vector<Real> y0{0.0};       // consistent because x0 + y0 = 1

DAESolverConfig config;
config.step_size = 0.05;
config.constraint_tol = 1e-10;

auto result = SolveDAEBackwardEuler(system, 0.0, x0, y0, 1.0, config);
auto xEnd = result.solution.getXValuesAtEnd();
auto yEnd = result.solution.getYValuesAtEnd();
```

The runnable demo obtains `x(1) = 0.574321814`, `y(1) = 0.425678186`, and final constraint norm
zero. Always inspect `result.status`, `failure_reason`, and constraint diagnostics before using a
trajectory.

### 25.2 Compute consistent initial algebraic variables

**Problem:** the differential initial state is known, but the algebraic variables are only a guess.
Use the consistent-IC helper to project the algebraic variables onto the constraint manifold.

```cpp
LinearDAE system;
Vector<Real> x0{0.8};
Vector<Real> yGuess{0.0};

auto ic = ComputeConsistentICDetailed(system, 0.0, x0, yGuess, 20, 1e-12);
bool verified = VerifyConsistentIC(system, 0.0, x0, yGuess, 1e-10);
```

For $x_0=0.8$, the helper updates `yGuess` to `0.2`, gives zero residual, and verification passes.
This is the right front door for measured or manually prepared DAE initial states: solve constraints
first, then integrate.

### 25.3 Compare Backward Euler, BDF2/BDF4, and Radau IIA

**Problem:** compare the built-in implicit methods on the same system and same step size.

```cpp
DAESolverConfig config;
config.step_size = 0.05;
config.constraint_tol = 1e-10;

auto backwardEuler = SolveDAEBackwardEuler(system, 0.0, x0, y0, 1.0, config);
auto bdf2 = SolveDAEBDF2(system, 0.0, x0, y0, 1.0, config);
auto bdf4 = SolveDAEBDF4(system, 0.0, x0, y0, 1.0, config);
auto radau = SolveDAERadauIIA(system, 0.0, x0, y0, 1.0, config);
```

At this step size, the final `x` values are about `0.574321814` (Backward Euler), `0.5677280454`
(BDF2), `0.5677036635` (BDF4), and `0.5676676418` (Radau IIA). The exact value is
$0.5+0.5e^{-2}\approx0.5676676416$, so the higher-order methods line up quickly. Constraint
preservation remains a separate diagnostic: the demo reports Radau's maximum constraint violation as
zero.

Use Backward Euler as a sturdy baseline, BDF methods for efficient stiff trajectories, and Radau IIA
when high-order stiff accuracy is worth the heavier coupled implicit solve.

### 25.4 Handle DAE events and restart with consistent constraints

**Problem:** stop at a zero-crossing event, mutate the state, and continue only after restoring the
algebraic constraint.

```cpp
class TimedDAEEvent : public IODESystemDAEWithEvents {
  Real eventTime;
public:
  explicit TimedDAEEvent(Real t) : eventTime(t) {}

  int getNumEvents() const override { return 1; }
  Real eventFunction(int, Real t, const Vector<Real>&, const Vector<Real>&) const override {
    return t - eventTime;
  }
  EventDirection getEventDirection(int) const override { return EventDirection::Increasing; }
  EventAction getEventAction(int) const override { return EventAction::Restart; }

  void handleEvent(int, Real, Vector<Real>& x, Vector<Real>& y) const override {
    x[0] = 0.25;
    y[0] = 42.0; // deliberately inconsistent; the solver recomputes y
  }
};
```

```cpp
TimedDAEEvent system(0.5);
Vector<Real> x0{0.0};
Vector<Real> y0{1.0};

DAESolverConfig config;
config.step_size = 0.2;
config.constraint_tol = 1e-10;

auto result = SolveDAEBackwardEulerWithEvents(system, 0.0, x0, y0, 1.0, config);
```

The demo detects one event at `t = 0.5`. `handleEvent` deliberately leaves `y` inconsistent, and the
event solver recomputes a consistent algebraic state before continuing; the final state is `x = 0.75`,
`y = 0.25`, with zero final constraint norm.

**Notes**

- MML's direct DAE solvers target semi-explicit index-1 systems. Higher-index problems need reformulation
  or index reduction before using this interface.
- Analytic Jacobians are not just performance polish; they make Newton solves and singular-Jacobian
  diagnostics much more reliable.
- A successful time integration is not enough by itself. Track `status`, Newton iteration counts,
  residual norms, and `max_constraint_violation`.
- Runnable version: `Cookbook_Recipe25_DAESolvers()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 26: Algebra: finite structures, symmetry, and exact arithmetic

MML's algebra layer lives in `<mml/base/Algebra_base.h>`. It covers finite permutations,
finite groups and group actions, modular rings, prime fields, exact polynomial arithmetic,
finite-field matrices, and small extension fields.

These examples deliberately focus on finite and exact algebra. Continuous rotation and rigid-motion
groups also live under Algebra, but rotations are introduced in [Recipe 17](#recipe-17-quaternions)
and fit more naturally with geometry and coordinate-transformation recipes.

### 26.1 Reorder real data with permutations

**Problem:** four hardware channels are wired in the wrong order. Use a permutation to move the
labels into their corrected positions, then use the inverse to map back to the raw layout.

**Setup**

```cpp
#include <mml/base/Algebra_base.h>

std::array<std::string, 4> wiredChannels = {
  "pressure", "temperature", "humidity", "voltage"
};
Algebra::Permutation<4> wiringFix = Algebra::Permutation<4>::from_cycles({ {0, 2, 1} });
```

**Solve**

```cpp
auto correctedChannels = wiringFix.apply(wiredChannels);
auto rawChannelsAgain = wiringFix.inverse().apply(correctedChannels);

auto images = wiringFix.images();
std::size_t order = wiringFix.order();
int sign = wiringFix.sign();
```

`images[i]` is the destination of index `i`. Applying the inverse recovers the original channel
layout. `order()` tells how many repeated applications return to the identity.

### 26.2 Model polygon symmetries with a dihedral group

**Problem:** a square marker can be rotated and mirrored. Work with the full symmetry group instead
of hand-writing separate rotation/reflection cases.

**Setup**

```cpp
Algebra::DihedralGroup square(4);
auto rotation = square.rotation(1);
auto reflection = square.reflection();
```

**Solve**

```cpp
auto relation = square.compose(square.compose(reflection, rotation), reflection);
auto rotationInverse = square.inverse(rotation);
auto reflectedVertices = square.permutation(reflection).images();
auto generated = Algebra::GeneratedSubgroup(square, { rotation, reflection });
```

The relation `s*r*s == r^-1` is the defining interaction between reflection and rotation in a
dihedral group. `permutation(element)` shows the action on polygon vertices, while
`GeneratedSubgroup` confirms that one rotation and one reflection generate all of `D4`.

### 26.3 Count distinct bracelets with Burnside's lemma

**Problem:** count visually distinct binary bracelets with six beads. There are $2^6=64$ raw
colorings, but many are the same after rotation or reflection.

**Setup**

```cpp
Algebra::DihedralGroup braceletSymmetry(6);
std::vector<std::vector<int>> colorings;
for (int mask = 0; mask < 64; ++mask) {
  std::vector<int> coloring(6);
  for (int position = 0; position < 6; ++position)
    coloring[position] = (mask >> position) & 1;
  colorings.push_back(coloring);
}

auto action = Algebra::MakeGroupAction<Algebra::DihedralElement, std::vector<int>>(
  [&braceletSymmetry](const Algebra::DihedralElement& element,
            const std::vector<int>& coloring) {
    return braceletSymmetry.permutation(element).apply(coloring);
  });
```

**Solve**

```cpp
std::size_t distinctBracelets = Algebra::BurnsideCount(braceletSymmetry, action, colorings);
std::size_t fixedByOneStepRotation = Algebra::FixedPointCount(
  braceletSymmetry, action, braceletSymmetry.rotation(1), colorings);
auto orbits = Algebra::OrbitPartition(braceletSymmetry, action, colorings);
```

Burnside's lemma computes

$$
|X/G|=\frac{1}{|G|}\sum_{g\in G}|\operatorname{Fix}(g)|.
$$

`BurnsideCount` also verifies that every transformed coloring remains in the supplied finite set,
which catches incomplete object lists.

### 26.4 Work exactly with modular arithmetic and prime fields

**Problem:** compare clock arithmetic with arithmetic in a prime field. Both are modular, but only
the prime field guarantees inverses for every nonzero element.

**Setup**

```cpp
using Z12 = Algebra::ModInt<12>;
using F7 = Algebra::PrimeFieldElement<7>;
Algebra::ModularRing<12> clockRing;
Algebra::PrimeField<7> field7;
```

**Solve**

```cpp
Z12 clockValue = Z12(10) + Z12(5);   // 3 mod 12
F7 quotient = F7(3) / F7(2);         // exact division in F7
F7 primitive = Algebra::PrimitiveRoot<7>();

bool ringLaws = Algebra::CheckRingLaws(clockRing);
bool fieldLaws = Algebra::CheckFieldLaws(field7);
```

`Z12(6).inverse()` throws because 6 is not a unit modulo 12. In `F7`, every nonzero element is a
unit, so division is exact finite-field arithmetic rather than floating-point arithmetic.

### 26.5 Solve an exact linear system over a finite field

**Problem:** solve a tiny modular constraint system over $F_7$ exactly, with no rounding or
conditioning questions.

**Setup**

```cpp
using F7 = Algebra::PrimeFieldElement<7>;
MatrixNM<F7, 2, 2> system{{ F7(2), F7(3) }, { F7(4), F7(1) }};
std::array<F7, 2> rightHandSide = { F7(1), F7(6) };
```

**Solve**

```cpp
auto inverse = Algebra::FieldMatrixInverse(system);
F7 determinant = Algebra::FieldMatrixDeterminant(system);

std::array<F7, 2> solution = {
  inverse(0, 0) * rightHandSide[0] + inverse(0, 1) * rightHandSide[1],
  inverse(1, 0) * rightHandSide[0] + inverse(1, 1) * rightHandSide[1]
};
```

`FieldMatrixInverse` performs Gaussian elimination over the field type itself. This is the right
tool for exact finite-field systems; do not use floating-point matrix inversion for this job.

### 26.6 Build a tiny extension field from polynomials

**Problem:** construct $GF(4)$ as $F_2[x]/(x^2+x+1)$ and verify the defining relation for
`alpha = x`.

**Setup**

```cpp
using F2 = Algebra::PrimeFieldElement<2>;

struct GF4Modulus {
  static Algebra::Polynomial<F2> modulus()
  {
    return Algebra::Polynomial<F2>({ F2(1), F2(1), F2(1) });
  }
};

using GF4 = Algebra::FiniteFieldElement<2, 2, GF4Modulus>;
GF4 alpha{ F2(0), F2(1) };
Algebra::ExtensionField<2, 2, GF4Modulus> field;
```

**Solve**

```cpp
GF4 alphaSquared = alpha * alpha;
GF4 alphaPlusOne = alpha + GF4(1);
GF4 inverse = alpha.inverse();

bool irreducible = Algebra::IsIrreducible(GF4Modulus::modulus());
bool fieldLaws = Algebra::CheckFieldLaws(field);
```

The modulus $x^2+x+1$ implies $x^2=x+1$ in characteristic 2, so `alpha * alpha` equals
`alpha + 1`. This is the same construction pattern used for larger finite fields, just small
enough to inspect directly.

**Notes**

- Permutations use zero-based images: `images[i]` is the destination of index `i`.
- `compose(left,right)` means `left * right`; for transformations, `right` is applied first.
- Finite law checks enumerate all elements, so they are meant for small structures and examples.
- `ModInt<n>` is a modular ring element; `PrimeFieldElement<p>` is the field-safe alias when `p`
  is prime.
- Extension-field modulus polynomials must be monic, irreducible, and have the declared degree.
- Runnable version: `Cookbook_Recipe26_Algebra()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 27: Complex analysis

The core API is `<mml/core/ComplexAnalysis.h>`, using `ComplexFunctionFromStdFunc` for
$f:\mathbb C\to\mathbb C$. Contours are parameterized as $\gamma(t)$ and integrated numerically
as $\int f(\gamma(t))\gamma'(t)\,dt$ with adaptive complex Simpson quadrature.

### 27.1 Differentiate an analytic complex function

**Problem:** differentiate $f(z)=e^z\sin z$ at $z=0.7+0.4i$ and compare finite-difference orders.

```cpp
#include <mml/core/ComplexAnalysis.h>

ComplexFunctionFromStdFunc f([](Complex z) { return std::exp(z) * std::sin(z); });
Complex z(0.7, 0.4);
Complex exact = std::exp(z) * (std::sin(z) + std::cos(z));

Complex secondOrder = Derivation::NDer2Complex(f, z);
Complex fourthOrder = Derivation::NDer4Complex(f, z);

DerivativeConfig config;
config.estimate_error = true;
auto detailed = Derivation::NDer4ComplexDetailed(f, z, config);
```

Errors are about $6.1\times10^{-11}$ and $7.4\times10^{-13}$, with a fourth-order estimate of
$1.3\times10^{-12}$. These are numerical derivatives in a real step direction; analyticity makes
the complex derivative direction-independent.

### 27.2 Integrate along complex contours

**Problem:** exercise every contour form and verify known antiderivatives and Cauchy's theorem.

```cpp
ComplexFunctionFromStdFunc identity([](Complex z) { return z; });
ComplexFunctionFromStdFunc one([](Complex) { return Complex(1,0); });
ComplexFunctionFromStdFunc square([](Complex z) { return z*z; });
ComplexFunctionFromStdFunc reciprocal([](Complex z) { return Complex(1,0) / z; });

ComplexAnalysis::LineSegmentContour line(Complex(0,0), Complex(1,1));
ComplexAnalysis::ArcContour arc(Complex(0,0), 1.0, 0.0, Constants::PI/2);
ComplexAnalysis::CircleContour circle(Complex(0,0), 1.0);

auto lineValue = ComplexAnalysis::ContourIntegral(identity, line);    // i
auto arcValue = ComplexAnalysis::ContourIntegral(one, arc);           // -1+i
auto analytic = ComplexAnalysis::ContourIntegral(square, circle);     // 0
auto singular = ComplexAnalysis::ContourIntegral(reciprocal, circle); // 2*pi*i
```

A custom parameterization supplies both $\gamma$ and $\gamma'$:

```cpp
auto ellipse = [](Real t) { return Complex(2*std::cos(t), std::sin(t)); };
auto ellipseDerivative = [](Real t) { return Complex(-2*std::sin(t), std::cos(t)); };
auto custom = ComplexAnalysis::ContourIntegral(
  identity, ellipse, ellipseDerivative, 0.0, 2*Constants::PI);
```

The custom closed integral is zero to about $4.4\times10^{-16}$. Circle contours are
counterclockwise; reversing parameter limits reverses the integral sign.

### 27.3 Compute winding numbers

**Problem:** determine whether a shifted circle surrounds selected points.

```cpp
ComplexAnalysis::CircleContour shifted(Complex(2,1), 3.0);
int aroundCenter = ComplexAnalysis::WindingNumber(shifted, Complex(2,1)); // 1
int aroundInside = ComplexAnalysis::WindingNumber(shifted, Complex(0,1)); // 1
int aroundOutside = ComplexAnalysis::WindingNumber(shifted, Complex(6,1)); // 0
```

Numerically this evaluates $(2\pi i)^{-1}\oint dz/(z-z_0)$ and rounds to an integer. The queried
point must not lie on the contour. MML's winding-number convenience API currently accepts circles.

### 27.4 Recover values and derivatives with Cauchy's formula

**Problem:** reconstruct $e^{z_0}$ and its first two derivatives using only values on a surrounding
circle.

```cpp
ComplexFunctionFromStdFunc exponential([](Complex z) { return std::exp(z); });
Complex z0(0.2, 0.1);
ComplexAnalysis::CircleContour contour(z0, 1.0);

Complex value = ComplexAnalysis::CauchyIntegralFormula(exponential, contour, z0);
Complex first = ComplexAnalysis::CauchyDerivative(exponential, contour, z0, 1);
Complex second = ComplexAnalysis::CauchyDerivative(exponential, contour, z0, 2);
```

All three equal $e^{z_0}$ to roughly $2.3\times10^{-16}$. The function must be analytic inside and
on the contour, and $z_0$ must be enclosed.

### 27.5 Compute residues and verify the residue theorem

**Problem:** analyze $f(z)=1/(z^2+1)$, whose simple poles are $\pm i$ with residues
$-i/2$ and $+i/2$.

```cpp
ComplexFunctionFromStdFunc f([](Complex z) {
  return Complex(1,0) / (z*z + Complex(1,0));
});

Complex residueI = ComplexAnalysis::Residue(f, Complex(0,1), 0.4);
Complex directI = ComplexAnalysis::ResidueSimplePole(f, Complex(0,1));
Complex residueMinusI = ComplexAnalysis::Residue(f, Complex(0,-1), 0.4);

auto onePole = ComplexAnalysis::ContourIntegral(
  f, ComplexAnalysis::CircleContour(Complex(0,1), 0.5));
auto bothPoles = ComplexAnalysis::ContourIntegral(
  f, ComplexAnalysis::CircleContour(Complex(0,0), 2.0));
```

The one-pole integral is $\pi$ because $2\pi i(-i/2)=\pi$. Enclosing both poles gives zero because
their residues cancel. A residue contour must isolate its singularity; `ResidueSimplePole` applies
only to first-order poles.

### 27.6 Count zeros and poles with the argument principle

**Problem:** count roots without first locating them, including multiplicity.

```cpp
ComplexFunctionFromStdFunc polynomial([](Complex z) { return z*z + Complex(1,0); });
int bothRoots = ComplexAnalysis::CountZeros(
  polynomial, ComplexAnalysis::CircleContour(Complex(0,0), 2.0)); // 2
int upperRoot = ComplexAnalysis::CountZeros(
  polynomial, ComplexAnalysis::CircleContour(Complex(0,1), 0.5)); // 1

ComplexFunctionFromStdFunc meromorphic([](Complex z) {
  return (z*z*z - Complex(1,0)) / (z - Complex(0.5,0));
});
int difference = ComplexAnalysis::ArgumentPrinciple(
  meromorphic, ComplexAnalysis::CircleContour(Complex(0,0), 2.0)); // 3-1 = 2
```

`CountZeros` assumes no poles inside. `ArgumentPrinciple` returns $N-P$, counting multiplicity.
The function must be meromorphic inside, with no zero or pole on the contour because $f'/f$ would
be singular there.

### 27.7 Count roots first, then locate them

**Problem:** establish that $z^3-1$ has three roots in a region, then numerically recover all three.

```cpp
ComplexFunctionFromStdFunc cubic([](Complex z) { return z*z*z - Complex(1,0); });
ComplexAnalysis::CircleContour contour(Complex(0,0), 2.0);
int expectedRoots = ComplexAnalysis::CountZeros(cubic, contour); // 3

auto root1 = RootFinding::FindRootNewtonComplex(cubic, Complex(1.1, 0.1));
auto root2 = RootFinding::FindRootNewtonComplex(cubic, Complex(-0.4, 0.8));
auto root3 = RootFinding::FindRootMuller(cubic, Complex(-0.4, -0.8));
```

The located roots are $1$, $-1/2+i\sqrt3/2$, and $-1/2-i\sqrt3/2$, with maximum residual about
$1.3\times10^{-12}$. Counting first tells you when root searches have missed a solution. Newton is
locally fast but guess-sensitive; Muller can naturally move through complex iterates.

**Notes**

- `ContourIntegrationResult` reports value, estimated error, evaluations, and convergence.
- Do not place singularities, zeros used by $f'/f$, or evaluation points on a contour.
- Smaller residue contours are not automatically better: extremely small radii amplify roundoff.
- Standard-library complex functions use principal branches. MML has no explicit branch-cut or
  analytic-continuation framework, so users must choose contours that respect the intended branch.
- Runnable version: `Cookbook_Recipe27_ComplexAnalysis()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).

---

## Recipe 30: Combinatorics and number theory

Exact counting and enumeration live in `<mml/base/Combinatorics.h>`; runtime 64-bit integer
number theory lives in `<mml/base/NumberTheory.h>`. Exact combinatorial values use `long long`
and throw `DomainError` on overflow. Counting may be fast even when enumerating every object is not.

### 30.1 Select a committee and assign roles

**Problem:** choose a four-person committee from ten people, then assign three distinct roles
(chair, secretary, treasurer) among its members.

```cpp
#include <mml/base/Combinatorics.h>

const int people = 10;
const int committeeSize = 4;
long long committees = Combinatorics::BinomialCoefficient(people, committeeSize); // 210
long long rolesPerCommittee = Combinatorics::FallingFactorial(committeeSize, 3);   // 24
long long totalAssignments = committees * rolesPerCommittee;                      // 5040
```

Enumeration uses zero-based person indices and callback visitors:

```cpp
long long enumerated = 0;
Combinatorics::ForEachCombination(people, committeeSize, [&](const std::vector<int>& committee) {
    Combinatorics::ForEachPermutation(committeeSize, [&](const std::vector<int>& roleOrder) {
        ++enumerated; // roleOrder[0..2] select committee positions for the three roles
    });
});
// enumerated = 5040
```

Use counting before enumeration: it predicts whether exhaustive generation is affordable.

### 30.2 Enumerate constrained resource allocations

**Problem:** distribute 12 identical resource units among exactly three projects, where every
project receives between 2 and 6 units. Project order does not matter.

```cpp
const int units = 12;
long long allPartitions = Combinatorics::PartitionCount(units); // 77
long long feasible = 0;

Combinatorics::ForEachPartition(units, [&](const std::vector<int>& allocation) {
    // Partitions arrive in non-increasing order.
    if (allocation.size() == 3 && allocation.front() <= 6 && allocation.back() >= 2)
        ++feasible;
});
// feasible = 5; first feasible partition is (6,4,2)
```

Integer partitions model allocations of indistinguishable units to unlabeled recipients. If project
identity matters, each distinct ordering must also be counted or generated.

### 30.3 Count ways to group labeled objects

**Problem:** split eight distinct students into nonempty teams.

```cpp
const int students = 8;
long long exactlyThreeUnlabeled = Combinatorics::StirlingSecond(students, 3); // 966
long long exactlyThreeNamed = exactlyThreeUnlabeled
    * Combinatorics::FallingFactorial(3, 3);                                  // 5796
long long anyNumberOfTeams = Combinatorics::BellNumber(students);             // 4140
```

$S(n,k)$ counts partitions of $n$ labeled objects into exactly $k$ nonempty **unlabeled** blocks.
Bell numbers sum over every possible nonempty block count. Multiplying by $k!$ labels the blocks.

### 30.4 Analyze an integer through its prime structure

**Problem:** factor an integer and derive arithmetic properties from its prime exponents.

```cpp
#include <mml/base/NumberTheory.h>

const long long value = 360360;
auto powers = NumberTheory::FactorizePowers(value);
// 2^3 * 3^2 * 5 * 7 * 11 * 13

long long divisors   = NumberTheory::DivisorCount(value);  // 192
long long divisorSum = NumberTheory::DivisorSum(value);    // 1572480
long long coprimes   = NumberTheory::EulerTotient(value);  // 69120
int mobius           = NumberTheory::MobiusMu(value);      // 0 (not square-free)
```

The same factorizer handles large 64-bit composites:

```cpp
std::uint64_t n = 1000000007ull * 1000000009ull;
bool prime = NumberTheory::IsPrime(n);       // false
auto factors = NumberTheory::Factorize(n);  // the two original primes
```

`IsPrime` uses deterministic Miller–Rabin over the full unsigned 64-bit range; factorization uses
Pollard rho and returns sorted factors with multiplicity.

### 30.5 Educational RSA-style encryption

**Problem:** connect primes, Euler's totient, modular inverse, and fast modular exponentiation in
a complete encryption/decryption round trip.

```cpp
const std::uint64_t p = 61, q = 53;
const std::uint64_t n = p * q;                         // 3233
const long long phi = (p - 1) * (q - 1);              // 3120
const std::uint64_t e = 17;
const std::uint64_t d = NumberTheory::ModInverse(e, phi); // 2753

const std::uint64_t message = 65;
std::uint64_t encrypted = NumberTheory::ModPow(message, e, n); // 2790
std::uint64_t decrypted = NumberTheory::ModPow(encrypted, d, n); // 65
```

This is an educational arithmetic demonstration, **not production cryptography**. Real RSA needs
large securely generated primes, padding, side-channel resistance, key validation, and audited code.

### 30.6 Synchronize repeating schedules with CRT

**Problem:** find the first time satisfying three repeating phase constraints:
$t\equiv2\pmod3$, $t\equiv3\pmod5$, and $t\equiv2\pmod7$.

```cpp
auto alignment = NumberTheory::ChineseRemainder({2, 3, 2}, {3, 5, 7});
// t = 23 mod 105: first alignment at 23, then every 105 time units
```

Non-coprime periods are also supported when their constraints agree modulo their gcd:

```cpp
auto sharedCycles = NumberTheory::ChineseRemainder({1, 4}, {6, 9});
// t = 13 mod 18

long long x, y;
long long g = NumberTheory::ExtendedGcd(6, 9, x, y);
// g = 3 = 6*(-1) + 9*(1); residue difference 4-1 is divisible by 3
```

An inconsistent system such as $t\equiv0\pmod2$, $t\equiv1\pmod4$ throws `DomainError`.
`ModInverse(a,m)` similarly requires $\gcd(a,m)=1$.

**Notes**

- `BinomialCoefficientReal` extends binomial range using floating-point log-gamma arithmetic; it is
  not an exact replacement for overflowing integer counts.
- `Fibonacci` is exact through $F_{92}$; other exact sequences also stop when 64-bit results overflow.
- Enumeration callbacks expose every object and therefore inherit the combinatorial explosion.
- `Factorize(0)` is outside the meaningful prime-factorization domain; positive arguments are the
  intended use for multiplicative functions.
- Runnable version: `Cookbook_Recipe30_CombinatoricsAndNumberTheory()` in
  [docs_demos_cookbook.cpp](../src/docs_demos/docs_demos_cookbook.cpp).
