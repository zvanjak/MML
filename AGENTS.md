# AI Agent Instructions for MinimalMathLibrary (MML)

> **Purpose:** This document provides AI assistants with comprehensive guidance for helping users leverage MML's numerical computing capabilities. Focus: **USAGE ONLY** — no contributing, no build systems, just pure mathematical power.

---

## 🎯 What is MML?

**MinimalMathLibrary (MML)** is a comprehensive, single-header C++ mathematical library for numerical computing. With just one `#include "MML.h"` directive, users get access to a complete toolkit — from vectors and matrices to ODE solvers and field operations.

**Key properties:**
- **Single-header library** — just `#include "MML.h"` and compute
- **Pure C++17** — no external dependencies
- **Cross-platform** — Windows, Linux, macOS
- **4,540 unit tests** — rigorously validated algorithms
- **Production-ready** — used in real scientific computing

---

## 📚 Core Capabilities Overview

| Domain | What MML Can Do |
|--------|-----------------|
| **Vectors & Matrices** | Full linear algebra, multiple coordinate systems, specialized storage |
| **Tensors** | Rank 2-5 tensors for advanced computations |
| **Linear Systems** | LU, QR, SVD, Cholesky solvers; handles ill-conditioned problems |
| **Eigenvalues** | Symmetric and general matrix eigensolvers |
| **Calculus** | Numerical derivatives (1st-8th order), integration (1D/2D/3D), improper integrals |
| **Root Finding** | Bisection, Newton, Brent, Secant methods |
| **ODE Solvers** | Fixed-step and adaptive solvers (RK4, Dormand-Prince, Bulirsch-Stoer) |
| **Optimization** | 1D (Golden Section, Brent), N-D (Nelder-Mead), Heuristic (SA, GA) |
| **Field Operations** | Gradient, divergence, curl, Laplacian in Cartesian/spherical/cylindrical |
| **Path & Surface Integrals** | Line integrals, surface integrals, flux calculations |
| **Differential Geometry** | Curvature, torsion, Frenet frames for parametric curves |
| **Dynamical Systems** | Fixed points, Lyapunov exponents, bifurcation analysis |
| **Interpolation** | Linear, polynomial, spline interpolation |
| **Statistics** | Descriptive stats, distributions, random number generation |
| **FFT/Signal Processing** | Forward/inverse FFT, power spectrum, real-valued transforms |
| **Geometry** | 2D/3D primitives, convex hull, triangulation, KD-trees |
| **Quaternions** | Full quaternion algebra, 3D rotations, SLERP |
| **Polynomials** | Evaluation, arithmetic, root finding, interpolation |

---

## 🚀 Getting Started Pattern

Always start with this minimal setup:

```cpp
#include "MML.h"
using namespace MML;

int main() {
    // MML code here
    return 0;
}
```

**Compile:** `g++ -std=c++17 -O3 myprogram.cpp -o myprogram`

---

## 📐 Vectors & Matrices

### Vector Creation
```cpp
// Dynamic vectors
Vector<Real> v(n);                        // Zero-initialized n-dimensional
Vector<Real> v{1.0, 2.0, 3.0};           // Initialize from list
Vector<Real> e = Vector<Real>::GetUnitVector(n, i);  // Unit vector e_i

// Fixed-size vectors (stack-allocated, faster)
VectorN<Real, 3> v{1, 2, 3};             // 3D vector
Vector2Cartesian v(x, y);                 // 2D Cartesian
Vector3Cartesian v(x, y, z);              // 3D Cartesian
Vector3Spherical v(r, theta, phi);        // 3D spherical (r, θ, φ)
Vector3Cylindrical v(rho, phi, z);        // 3D cylindrical (ρ, φ, z)
```

### Vector Operations
```cpp
v + w, v - w, v * scalar, v / scalar      // Arithmetic
v.NormL1(), v.NormL2(), v.NormLInf()      // Norms: L1, L2 (Euclidean), L∞

// Dot and cross products (3D Cartesian vectors)
Vector3Cartesian a{1, 2, 3}, b{4, 5, 6};
Real dot = ScalarProduct(a, b);           // a·b = 32
Vector3Cartesian cross = VectorProduct(a, b);  // a×b
```

### Matrix Creation
```cpp
Matrix<Real> A(m, n);                     // m×n zero matrix
Matrix<Real> A{3, 3, {1,2,3, 4,5,6, 7,8,9}};  // Row-major initialization
Matrix<Real> I = Matrix<Real>::Identity(n);   // n×n identity

// Specialized matrices
MatrixSym<Real> S(n);                     // Symmetric matrix (stores only upper triangle)
MatrixTriDiag<Real> T(n);                 // Tridiagonal matrix
MatrixNM<Real, 3, 3> M;                   // Fixed-size 3×3 (stack-allocated)
```

### Matrix Operations
```cpp
A + B, A - B, A * B                       // Matrix arithmetic
A * v                                     // Matrix-vector product
A.GetTranspose()                          // Transpose Aᵀ
A.GetInverse()                            // Inverse A⁻¹
A.Determinant()                           // Determinant |A|
A.Trace()                                 // Trace Σ Aᵢᵢ
A.NormFrobenius()                         // Frobenius norm
A.GetRow(i), A.GetColumn(j)               // Extract row/column
A.Submatrix(r1, r2, c1, c2)               // Extract submatrix
```

---

## ⚙️ Linear Systems (Ax = b)

### Solver Selection Guide

| Solver | Use When | Performance |
|--------|----------|-------------|
| **LUSolver** | General square matrices, multiple RHS | Fast, O(n³) |
| **CholeskySolver** | Symmetric positive definite (SPD) | Fastest, 2× faster than LU |
| **QRSolver** | Ill-conditioned, overdetermined systems | Stable, O(n³) |
| **SVDecompositionSolver** | Rank-deficient, singular, near-singular | Most robust |
| **GaussJordanSolver** | Small systems, need matrix inverse | Simple |

### Usage Examples

```cpp
Matrix<Real> A{3, 3, {4, 1, 2,  1, 5, 1,  2, 1, 6}};
Vector<Real> b{8, 9, 12};

// LU Decomposition (general purpose)
LUSolver<Real> lu(A);
Vector<Real> x = lu.Solve(b);
Real det = lu.det();                      // Determinant
Matrix<Real> Ainv; lu.inverse(Ainv);      // Matrix inverse

// QR Decomposition (overdetermined/least-squares)
QRSolver<Real> qr(A);
Vector<Real> x = qr.Solve(b);             // Square system
Vector<Real> x = qr.LeastSquaresSolve(b); // Minimize ||Ax-b||₂

// SVD (rank-deficient, singular)
SVDecompositionSolver svd(A);
Vector<Real> x = svd.Solve(b);
int rank = svd.Rank();                    // Matrix rank
Real cond = 1.0/svd.inv_condition();      // Condition number

// Cholesky (SPD matrices - fastest)
CholeskySolver<Real> chol(A);             // A must be SPD!
Vector<Real> x = chol.Solve(b);
```

### Eigenvalues & Eigenvectors

```cpp
// Symmetric matrices (real eigenvalues, orthogonal eigenvectors)
auto result = SymmMatEigenSolverJacobi::Solve(A);
Vector<Real> eigenvals = result.eigenvalues;
Matrix<Real> eigenvecs = result.eigenvectors;

// Access individual
Real lambda_i = result.eigenvalues[i];
Vector<Real> v_i = result.eigenvectors.VectorFromColumn(i);
```

---

## 📈 Calculus

### Numerical Derivatives

**Accuracy selection guide:** Higher order = more accuracy but more function evaluations.

| Method | Order | Accuracy | Function Evals | Use For |
|--------|-------|----------|----------------|---------|
| `NDer1` | 1 | O(h) | 2 | Quick estimates |
| `NDer2` | 2 | O(h²) | 3 | General use (recommended) |
| `NDer4` | 4 | O(h⁴) | 5 | High accuracy |
| `NDer6` | 6 | O(h⁶) | 7 | Very high accuracy |
| `NDer8` | 8 | O(h⁸) | 9 | Near machine precision |

```cpp
RealFunction f = [](Real x) { return std::sin(x); };

// First derivative f'(x₀)
Real df = Derivation::NDer2(f, x0);       // Recommended for most cases
Real df = Derivation::NDer4(f, x0);       // High accuracy
Real df = Derivation::NDer8(f, x0);       // Maximum precision

// Second derivative f''(x₀)
Real d2f = Derivation::NSecDer4(f, x0);

// Third derivative f'''(x₀)
Real d3f = Derivation::NThirdDer4(f, x0);

// With error estimate
Real error;
Real df = Derivation::NDer2(f, x0, &error);

// Gradient of multivariate function f: ℝⁿ → ℝ
ScalarFunction<3> g([](const VectorN<Real, 3>& v) {
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
});
auto gradient = Derivation::Gradient<3>(g, point);  // Returns VectorN<Real, 3>
```

### Numerical Integration (1D)

**Method selection guide:**

| Method | Use For | Accuracy |
|--------|---------|----------|
| `IntegrateRomberg` | Smooth functions | High (Richardson extrapolation) |
| `IntegrateSimpson` | General purpose | Good |
| `IntegrateTrap` | Quick estimates | Basic |
| `IntegrateGaussKronrod` | High precision needed | Very high |
| `GaussLegendre` | Polynomial-like functions | Excellent for smooth |
| `GaussLaguerre` | e⁻ˣ weight, [0,∞) | Specialized |
| `GaussHermite` | e⁻ˣ² weight, (-∞,∞) | Specialized |

```cpp
RealFunction f = [](Real x) { return std::sin(x); };

// Adaptive quadrature (recommended)
auto result = IntegrateRomberg(f, a, b);  // Returns IntegrationResult
Real I = result.value;
bool converged = result.converged;
Real error = result.error_estimate;

// Gauss-Kronrod (highest accuracy)
auto result = IntegrateGaussKronrod(f, a, b, GK_G10K21);

// Improper integrals
auto result = IntegrateUpperInf(f, a);         // ∫ₐ^∞ f(x)dx
auto result = IntegrateLowerInf(f, b);         // ∫₋∞^b f(x)dx
auto result = IntegrateInfInf(f);              // ∫₋∞^∞ f(x)dx
```

### Multi-dimensional Integration

```cpp
// 2D integration over rectangle
ScalarFunction<2> f2d([](const VectorN<Real, 2>& x) {
    return x[0]*x[0] + x[1]*x[1];
});
auto result = Integrate2D(f2d, GAUSS10, ax, bx, ay, by, nx, ny);

// 3D integration
ScalarFunction<3> f3d([](const VectorN<Real, 3>& x) { ... });
auto result = Integrate3D(f3d, GAUSS10, ax, bx, y_lo, y_hi, z_lo, z_hi);

// Variable limits (e.g., triangular region)
auto y_lo = [](Real x) { return 0.0; };
auto y_hi = [](Real x) { return 1.0 - x; };  // y from 0 to 1-x
auto z_lo = [](Real x, Real y) { return 0.0; };
auto z_hi = [](Real x, Real y) { return 1.0 - x - y; };

// Monte Carlo (for high dimensions)
MonteCarloIntegrator<N> mc;
auto result = mc.integrate(func, lower, upper, config.samples(100000));
```

---

## 🎯 Root Finding

### 1D Root Finding

| Method | Requires | Convergence | Best For |
|--------|----------|-------------|----------|
| `FindRootBisection` | Bracket [a,b] | Linear | Guaranteed convergence |
| `FindRootBrent` | Bracket [a,b] | Superlinear | **Recommended general use** |
| `FindRootNewton` | Initial guess | Quadratic | Fast when derivative easy |
| `FindRootSecant` | Two guesses | Superlinear | No derivative available |
| `FindRootRidders` | Bracket [a,b] | Quadratic | Alternative to Brent |

```cpp
RealFunction f = [](Real x) { return x*x - 2; };  // Root: x = √2

// Brent's method (recommended)
Real root = RootFinding::FindRootBrent(f, 1.0, 2.0, 1e-12);

// Newton-Raphson (faster if derivative available)
Real root = RootFinding::FindRootNewton(f, 1.0, 2.0, 1e-12);

// Find multiple roots via bracketing
Vector<Real> lo, hi;
int n = RootFinding::FindRootBrackets(f, -10.0, 10.0, 100, lo, hi);
```

### Polynomial Roots

```cpp
// Quadratic: ax² + bx + c = 0
Complex r1, r2;
int numReal = SolveQuadratic(a, b, c, r1, r2);

// Cubic: ax³ + bx² + cx + d = 0
Complex c1, c2, c3;
SolveCubic(a, b, c, d, c1, c2, c3);

// Quartic: ax⁴ + bx³ + cx² + dx + e = 0
Complex q1, q2, q3, q4;
SolveQuartic(a, b, c, d, e, q1, q2, q3, q4);
```

---

## 🌊 Differential Equations (ODEs)

### Defining ODE Systems

```cpp
// Method 1: Lambda-based (simplest)
ODESystem system([](Real t, const Vector<Real>& y, Vector<Real>& dydt) {
    dydt[0] = y[1];           // dx/dt = v
    dydt[1] = -y[0];          // dv/dt = -x (harmonic oscillator)
}, 2);  // dimension = 2

// Method 2: Class-based (for complex systems)
class MyODE : public IODESystem {
    int getDim() const override { return 2; }
    void derivs(Real t, const Vector<Real>& y, Vector<Real>& dydt) const override {
        dydt[0] = y[1];
        dydt[1] = -y[0];
    }
};
```

### Solver Selection Guide

| Solver | Type | Use For |
|--------|------|---------|
| `ODESystemStepperEuler` | Fixed-step | Educational only |
| `ODESystemStepperRK4` | Fixed-step | General problems, simple |
| `DormandPrince5_Stepper` | Adaptive | **Recommended general use** |
| `DormandPrince8_Stepper` | Adaptive | High accuracy requirements |
| `CashKarp_Stepper` | Adaptive | Alternative to DP5 |
| `BulirschStoer_Stepper` | Adaptive | Very smooth solutions |

```cpp
Vector<Real> y0{1.0, 0.0};  // Initial conditions

// Fixed-step solver
ODESystemStepperRK4 stepper;
ODESystemFixedStepSolver solver(system, stepper);
ODESystemSolution sol = solver.integrate(y0, t0, tend, numSteps);

// Adaptive solver (recommended)
ODEAdaptiveIntegrator<DormandPrince5_Stepper> solver(system);
ODESystemSolution sol = solver.integrate(y0, t0, tend, outputInterval, eps);

// Access solution
Real t_i = sol.t(i);              // Time at step i
Real y_j_i = sol.y(i, j);         // Component j at step i
Vector<Real> y_i = sol.yVec(i);   // Full state at step i
```

### Predefined Dynamical Systems

MML includes famous dynamical systems ready to use:

```cpp
LorenzSystem lorenz(sigma, rho, beta);    // Lorenz attractor
RosslerSystem rossler(a, b, c);           // Rössler system
VanDerPolSystem vdp(mu);                  // Van der Pol oscillator
DoublePendulumSystem dp(L1, L2, g);       // Double pendulum
DuffingSystem duffing(alpha, beta, delta);// Duffing oscillator
HodgkinHuxleySystem hh(...);              // Neuron model
```

---

## 📉 Optimization

### 1D Minimization

```cpp
RealFunction f = [](Real x) { return x*x + 3*x + 2; };

// Brent's method (recommended)
Real xmin = Optimization::Brent(f, a, b, tol);

// Golden section
Real xmin = Optimization::GoldenSection(f, a, b, tol);
```

### Multi-dimensional Minimization

```cpp
// Define objective function
class MyFunc : public IScalarFunction<N> {
    Real operator()(const VectorN<Real, N>& x) const override {
        return x[0]*x[0] + x[1]*x[1];  // Rosenbrock, sphere, etc.
    }
};

MyFunc func;
VectorN<Real, N> x0{...};  // Initial guess

// Nelder-Mead simplex (derivative-free)
NelderMead<N> optimizer(maxIter, tol);
auto result = optimizer.Minimize(func, x0, delta);

// Simulated annealing (global optimization)
SimulatedAnnealing<N> sa(config);
auto result = sa.Minimize(func, x0);

// Genetic algorithm (global optimization)
GeneticAlgorithm<N> ga(config);
auto result = ga.Minimize(func, x0);
```

---

## 🌀 Field Operations

### Scalar Fields (f: ℝⁿ → ℝ)

```cpp
ScalarFunction<3> phi([](const VectorN<Real, 3>& x) {
    return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];  // r²
});
VectorN<Real, 3> point{1, 2, 3};

// Gradient ∇f (vector pointing uphill)
auto grad = ScalarFieldOperations::GradientCart<3>(phi, point);

// Laplacian ∇²f (divergence of gradient)
Real lap = ScalarFieldOperations::LaplacianCart<3>(phi, point);

// In other coordinate systems
auto grad_sph = ScalarFieldOperations::GradientSpher(phi, point);
auto grad_cyl = ScalarFieldOperations::GradientCyl(phi, point);
```

### Vector Fields (F: ℝⁿ → ℝⁿ)

```cpp
VectorFunction<3> F([](const VectorN<Real, 3>& x) -> VectorN<Real, 3> {
    return {x[1], -x[0], x[2]};  // Rotating field
});

// Divergence ∇·F (expansion rate)
Real div = VectorFieldOperations::DivCart<3>(F, point);

// Curl ∇×F (rotation/vorticity) - 3D only
auto curl = VectorFieldOperations::CurlCart(F, point);

// Jacobian matrix
auto J = VectorFieldOperations::JacobianCart<3>(F, point);
```

---

## 🛤️ Path & Surface Integrals

### Line Integrals

```cpp
// Vector field
VectorFunction<3> F([](const VectorN<Real, 3>& p) -> VectorN<Real, 3> {
    return {p[1], -p[0], 0.0};
});

// Parametric curve
ParametricCurve<3> circle([](Real t) -> VectorN<Real, 3> {
    return {cos(t), sin(t), 0.0};
});

// Work integral: ∫ F·dr
Real work = PathIntegration::LineIntegral(F, circle, 0.0, 2*Constants::PI, tol);

// Scalar path integral: ∫ f ds (arc length weighted)
ScalarFunction<3> density([](const VectorN<Real, 3>& p) { return 1.0; });
Real arc_length = PathIntegration::ScalarPathIntegral(density, circle, 0, 2*PI, tol);
```

### Surface Integrals

```cpp
// Flux through surface: ∮∮ F·n̂ dS
Cube3D cube(1.0, Point3Cartesian(0.5, 0.5, 0.5));
Real flux = SurfaceIntegration::SurfaceIntegral(F, cube, tol);

// For parametric surfaces
ParametricSurface<3> sphere([](Real u, Real v) -> VectorN<Real, 3> {
    return {sin(u)*cos(v), sin(u)*sin(v), cos(u)};
});
Real flux = SurfaceIntegration::SurfaceIntegral(F, sphere, u_lo, u_hi, v_lo, v_hi, tol);
```

---

## 🔄 Coordinate Transformations

```cpp
// Convert between coordinate systems
Vector3Cartesian cart{1.0, 1.0, 1.0};
Vector3Spherical sph = CoordTransfCartToSpher.transf(cart);
Vector3Cylindrical cyl = CoordTransfCartToCyl.transf(cart);

// Convert back
Vector3Cartesian back = CoordTransfSpherToCart.transf(sph);

// Transform vectors (covariant for forces/gradients)
auto force_cart = CoordTransfSpherToCart.transfVecCovariant(force_sph, point);

// Transform vectors (contravariant for velocities)
auto vel_cart = CoordTransfSpherToCart.transfVecContravariant(vel_sph, point);
```

---

## 📊 Parametric Curves & Differential Geometry

```cpp
// Define 3D curve: helix r(t) = (cos(t), sin(t), 0.2t)
ParametricCurve<3> helix([](Real t) -> VectorN<Real, 3> {
    return {cos(t), sin(t), 0.2*t};
});

// Curve properties at parameter t
auto pos = helix(t);                    // Position r(t)
auto tangent = helix.getTangent(t);     // Tangent dr/dt
auto unit_tan = helix.getTangentUnit(t);// Unit tangent T
auto normal = helix.getNormal(t);       // Principal normal N
auto binormal = helix.getBinormal(t);   // Binormal B = T × N

// Curvature κ and torsion τ (Frenet-Serret)
Real kappa = helix.getCurvature(t);
Real tau = helix.getTorsion(t);

// Arc length from t=a to t=b
Real length = helix.getArcLength(a, b);

// Predefined curves
Curves::HelixCurve helix(radius, pitch);
Curves::LemniscateCurve lemniscate;
Curves::ToroidalSpiralCurve torus(R, r);
```

---

## 🦋 Dynamical Systems Analysis

### Fixed Point Analysis

```cpp
LorenzSystem lorenz(10.0, 28.0, 8.0/3.0);

// Find fixed points
std::vector<Vector<Real>> guesses = {
    Vector<Real>{0, 0, 0},
    Vector<Real>{8, 8, 27}
};
auto fixedPoints = FixedPointFinder::FindMultiple(lorenz, guesses);

for (const auto& fp : fixedPoints) {
    std::cout << "Location: " << fp.location << std::endl;
    std::cout << "Type: " << ToString(fp.type) << std::endl;  // Node, Spiral, Saddle
    std::cout << "Stable: " << fp.isStable << std::endl;
    // Eigenvalues determine stability
    for (const auto& ev : fp.eigenvalues)
        std::cout << "λ = " << ev.real() << " + " << ev.imag() << "i" << std::endl;
}
```

### Lyapunov Exponents

```cpp
// Quantify chaos: positive Lyapunov exponent = chaos!
auto result = LyapunovAnalyzer::Compute(
    lorenz,
    lorenz.getDefaultInitialCondition(),
    500.0,   // Integration time
    1.0,     // Orthonormalization interval
    0.01     // Step size
);

std::cout << "Lyapunov spectrum: " << result.exponents << std::endl;
std::cout << "Max exponent: " << result.maxExponent << std::endl;  // >0 = chaos
std::cout << "Kaplan-Yorke dimension: " << result.kaplanYorkeDimension << std::endl;
```

### Bifurcation Analysis

```cpp
// Sweep parameter to find bifurcations
auto bifurcation = BifurcationAnalyzer::Sweep(
    lorenz,
    1,              // Parameter index (rho)
    20.0, 30.0,     // Parameter range
    50,             // Number of steps
    Vector<Real>{1, 1, 1},  // Initial condition
    2,              // Record z-component maxima
    100.0, 50.0     // Transient, recording time
);
```

---

## 📈 Interpolation & Approximation

```cpp
Vector<Real> x_data{0, 1, 2, 3, 4};
Vector<Real> y_data{0, 1, 4, 9, 16};

// Create interpolating functions
LinearInterpRealFunc linear(x_data, y_data);    // Piecewise linear
SplineInterpRealFunc spline(x_data, y_data);    // Cubic spline (smooth)
PolynomInterpRealFunc poly(x_data, y_data, 3);  // Polynomial degree 3

// Evaluate at any point in range
Real y = spline(2.5);  // Interpolate at x=2.5
```

---

## 📊 Statistics & Probability

### Descriptive Statistics

```cpp
std::vector<Real> data{1.2, 3.4, 2.1, 5.6, ...};

Real mean = Statistics::Mean(data);
Real median = Statistics::Median(data);
Real variance = Statistics::Variance(data);
Real stddev = Statistics::StdDev(data);
Real skewness = Statistics::Skewness(data);
Real kurtosis = Statistics::Kurtosis(data);
```

### Probability Distributions

```cpp
// Normal distribution
Real pdf = NormalDistribution::PDF(x, mu, sigma);
Real cdf = NormalDistribution::CDF(x, mu, sigma);
Real quantile = NormalDistribution::InverseCDF(p, mu, sigma);

// Other distributions available:
// BinomialDistribution, PoissonDistribution, ExponentialDistribution,
// ChiSquareDistribution, StudentTDistribution, FDistribution, etc.
```

### Random Number Generation

```cpp
Ran rng(seed);                        // Initialize RNG

Real u = rng.doub();                  // Uniform [0,1)
int i = rng.int32(n);                 // Uniform [0, n)

NormalDeviate normal(rng, mu, sigma); // Normal distribution
Real x = normal.dev();

ExponentialDeviate exp(rng, lambda);  // Exponential distribution
PoissonDeviate poisson(rng, lambda);  // Poisson distribution
BinomialDeviate binom(rng, n, p);     // Binomial distribution
```

---

## 🔊 Fourier & Signal Processing

```cpp
Vector<Complex> signal(n);  // n must be power of 2

// Forward FFT
FFT::Transform(signal, 1);

// Inverse FFT (remember to normalize)
FFT::Transform(signal, -1);
for (auto& val : signal) val /= n;

// Real-valued FFT (2× faster for real signals)
Vector<Real> real_signal(n);
FFT::RealFFT(real_signal);    // Forward
FFT::RealIFFT(real_signal);   // Inverse

// Power spectrum
Vector<Real> power = FFT::PowerSpectrum(signal);
```

---

## 📐 2D & 3D Geometry

### 2D Primitives

```cpp
Point2Cartesian<Real> p(x, y);
Line2D<Real> line(p1, p2);
Circle2D<Real> circle(center, radius);
Polygon2D<Real> poly{{p1, p2, p3, ...}};

Real area = poly.Area();
bool inside = poly.Contains(point);
Real dist = line.Distance(point);
```

### 3D Primitives

```cpp
Point3Cartesian<Real> p(x, y, z);
Line3D<Real> line(p1, p2);
Plane3D<Real> plane(point, normal);
Triangle3D<Real> tri(p1, p2, p3);
Sphere3D<Real> sphere(center, radius);
Cube3D cube(size, center);

Real area = tri.Area();
Vector3Cartesian normal = tri.Normal();
Real dist = plane.Distance(point);
```

---

## 🧮 Polynomials

```cpp
// Create polynomial: p(x) = 2x³ - 3x² + x - 5
PolynomReal p{-5, 1, -3, 2};  // Coefficients: [constant, x, x², x³]

// Evaluate
Real val = p(2.0);

// Arithmetic
PolynomReal sum = p1 + p2;
PolynomReal prod = p1 * p2;

// Calculus
PolynomReal dp = p.derivative();      // p'(x)
PolynomReal ip = p.integral();        // ∫p(x)dx

// Create from data points (Lagrange interpolation)
Polynom<Real> p = Polynom<Real>::FromValues(x_vals, y_vals);
```

---

## 🔄 Quaternions (3D Rotations)

```cpp
Quaternion q(w, x, y, z);                    // Create quaternion
Quaternion q = Quaternion::FromAxisAngle(axis, angle);  // From rotation

// Operations
Quaternion product = q1 * q2;                 // Quaternion multiplication
Quaternion conjugate = q.Conjugate();
Quaternion inverse = q.Inverse();
Real norm = q.Norm();
Quaternion unit = q.Normalize();

// Rotate vector
Vector3Cartesian rotated = q.Rotate(v);

// Interpolation
Quaternion interp = Quaternion::Slerp(q1, q2, t);  // Spherical interpolation

// Convert to/from Euler angles
auto [roll, pitch, yaw] = q.ToEulerAngles();
Quaternion q = Quaternion::FromEulerAngles(roll, pitch, yaw);
```

---

## 💡 Best Practices for AI Agents

### Solver Recommendations

| Problem Type | Recommended Approach |
|--------------|---------------------|
| Linear system (general) | `LUSolver` |
| Linear system (SPD) | `CholeskySolver` (2× faster) |
| Linear system (ill-conditioned) | `QRSolver` or `SVDecompositionSolver` |
| Least squares | `QRSolver::LeastSquaresSolve` |
| Smooth integration | `IntegrateRomberg` |
| High-precision integration | `IntegrateGaussKronrod` with `GK_G10K21` |
| Root finding | `FindRootBrent` (no derivative) or `FindRootNewton` (with derivative) |
| ODE (general) | `ODEAdaptiveIntegrator<DormandPrince5_Stepper>` |
| ODE (high accuracy) | `DormandPrince8_Stepper` or `BulirschStoer_Stepper` |
| 1D optimization | `Optimization::Brent` |
| N-D optimization | `NelderMead` (local), `SimulatedAnnealing` (global) |

### Error Handling

MML uses exceptions for error conditions:

```cpp
try {
    LUSolver<Real> solver(singular_matrix);
    Vector<Real> x = solver.Solve(b);
} catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
}
```

### Precision Considerations

- Use `NDer4` or higher for derivatives of smooth functions
- Check `IntegrationResult::converged` for integration reliability
- For ill-conditioned matrices, check condition number via SVD
- ODE solutions should be validated against known analytical cases when possible

---

## 📚 Quick Reference Links

| Documentation | Contents |
|---------------|----------|
| `docs/API_CHEATSHEET.md` | Compact API reference |
| `docs/QUICK_START.md` | 5-minute getting started |
| `docs/base/` | Vectors, matrices, tensors, functions |
| `docs/core/` | Derivation, integration, linear solvers |
| `docs/algorithms/` | ODE, root finding, optimization |
| `docs/systems/` | Dynamical systems analysis |
| `examples/` | Working code examples |
| `tests/` | Unit tests (great usage examples) |

---

## 🎯 Common Task Patterns

### Task: Solve a physics simulation with ODEs

```cpp
// 1. Define the system
ODESystem system([](Real t, const Vector<Real>& y, Vector<Real>& dydt) {
    // Your physics equations here
    dydt[0] = y[1];
    dydt[1] = -9.81;  // Gravity
}, 2);

// 2. Set initial conditions
Vector<Real> y0{10.0, 0.0};  // Height=10, velocity=0

// 3. Integrate
ODEAdaptiveIntegrator<DormandPrince5_Stepper> integrator(system);
auto solution = integrator.integrate(y0, 0.0, 5.0, 0.01, 1e-8);

// 4. Analyze results
for (int i = 0; i < solution.size(); i++) {
    std::cout << "t=" << solution.t(i) << " y=" << solution.y(i, 0) << std::endl;
}
```

### Task: Find where two functions intersect

```cpp
RealFunction f = [](Real x) { return std::sin(x); };
RealFunction g = [](Real x) { return 0.5 * x; };
RealFunction diff = [&](Real x) { return f(x) - g(x); };

// Find roots of f(x) - g(x) = 0
Real intersection = RootFinding::FindRootBrent(diff, 0.5, 2.0, 1e-10);
```

### Task: Compute gradient and find minimum

```cpp
// Objective function
ScalarFunction<2> f([](const VectorN<Real, 2>& x) {
    return (x[0]-1)*(x[0]-1) + (x[1]-2)*(x[1]-2);  // Min at (1,2)
});

// Compute gradient at a point
VectorN<Real, 2> point{0.0, 0.0};
auto grad = Derivation::Gradient<2>(f, point);
std::cout << "Gradient: " << grad << std::endl;

// Find minimum
NelderMead<2> optimizer(1000, 1e-10);
auto result = optimizer.Minimize(f, point, 0.5);
std::cout << "Minimum at: " << result.min_point << std::endl;
```

### Task: Verify a vector calculus theorem

```cpp
// Gauss's Divergence Theorem: ∫∫∫(∇·F)dV = ∮∮(F·n̂)dS
VectorFunction<3> F([](const VectorN<Real, 3>& p) -> VectorN<Real, 3> {
    return {p[0]*p[0], p[1]*p[1], p[2]*p[2]};
});

// Divergence: ∇·F = 2x + 2y + 2z
ScalarFunction<3> divF([&F](const VectorN<Real, 3>& p) {
    return VectorFieldOperations::DivCart<3>(F, p);
});

// Volume integral over unit cube
auto y_lo = [](Real) { return 0.0; };  auto y_hi = [](Real) { return 1.0; };
auto z_lo = [](Real,Real) { return 0.0; };  auto z_hi = [](Real,Real) { return 1.0; };
Real vol = Integrate3D(divF, GAUSS10, 0, 1, y_lo, y_hi, z_lo, z_hi).value;

// Surface integral
Cube3D cube(1.0, Point3Cartesian(0.5, 0.5, 0.5));
Real surf = SurfaceIntegration::SurfaceIntegral(F, cube, 1e-8);

std::cout << "Volume integral:  " << vol << std::endl;
std::cout << "Surface integral: " << surf << std::endl;
std::cout << "Theorem verified: " << (std::abs(vol - surf) < 1e-10) << std::endl;
```

---

*MML v1.2 — Comprehensive C++ Numerical Computing Library*
*This document is for AI assistants helping users leverage MML's capabilities.*
