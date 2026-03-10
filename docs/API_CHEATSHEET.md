# MML API Quick Reference

**MinimalMathLibrary** — Fast reference for experienced users. Complete docs at [docs/](.).

```cpp
#include "MML.h"  // Single header includes everything
using namespace MML;
```

---

## 📊 Vectors

### Construction & Access
```cpp
Vector<Real> v(n);                    // Zero-initialized n-dimensional vector
Vector<Real> v{1.0, 2.0, 3.0};       // Initialize from list
Vector<Real> v = Vector<Real>::GetUnitVector(n, i);  // Unit vector e_i
v[i] = val;                           // Element access (0-indexed)
v.size();                             // Number of elements
```

### Operations
```cpp
v + w, v - w, v * scalar, v / scalar  // Arithmetic operations
v.NormL1();                           // L1 norm: Σ|vᵢ|
v.NormL2();                           // L2 norm (Euclidean): √(Σvᵢ²)
v.NormLInf();                         // L∞ norm: max|vᵢ|
```

**Note:** For dot product and cross product, use specialized vector types (`Vector3Cartesian`, etc.):
```cpp
Vector3Cartesian v{1, 2, 3}, w{4, 5, 6};
Real dot = ScalarProduct(v, w);       // Dot product: v·w
Vector3Cartesian cross = VectorProduct(v, w);  // Cross product: v×w
```

### Fixed-Size Vectors
```cpp
VectorN<Real, 3> v{1, 2, 3};         // Stack-allocated 3D vector
Vector2Cartesian v(x, y);             // 2D Cartesian
Vector3Cartesian v(x, y, z);          // 3D Cartesian
Vector3Spherical v(r, theta, phi);    // 3D spherical
Vector3Cylindrical v(rho, phi, z);    // 3D cylindrical
```

---

## 🔢 Matrices

### Construction & Access
```cpp
Matrix<Real> A(m, n);                 // m×n matrix (zero-initialized)
Matrix<Real> A{m, n, {1,2,3,4,...}};  // Initialize from list (row-major)
Matrix<Real> I = Matrix<Real>::Identity(n);  // n×n identity matrix
A(i, j) = val;                        // Element access (0-indexed)
A.RowNum(), A.ColNum();               // Dimensions
```

### Operations
```cpp
A + B, A - B, A * B                   // Matrix arithmetic
A * v                                 // Matrix-vector product
A * scalar                            // Scalar multiplication
A.transpose();                     // Transpose: Aᵀ
A.Trace();                            // Trace: Σ Aᵢᵢ
A.NormFrobenius();                    // Frobenius norm: √(Σ Aᵢⱼ²)
A.GetRow(i), A.GetColumn(j);          // Extract row/column vector
A.Submatrix(r1,r2, c1,c2);            // Extract submatrix
```

### Special Matrices
```cpp
MatrixNM<Real, 3, 3> M;               // Fixed-size 3×3 matrix (stack)
MatrixSym<Real> S(n);                 // Symmetric n×n matrix
MatrixTriDiag<Real> T(n);             // Tridiagonal matrix
```

---

## ⚙️ Linear Systems (A·x = b)

### Direct Solvers
```cpp
// LU Decomposition (general square systems, multiple RHS)
LUSolver<Real> lu(A);
Vector<Real> x = lu.Solve(b);
Real det = lu.det();                  // Determinant
Matrix<Real> Ainv; lu.inverse(Ainv);  // Matrix inverse

// QR Decomposition (ill-conditioned, overdetermined)
QRSolver<Real> qr(A);
Vector<Real> x = qr.Solve(b);         // Square system
Vector<Real> x = qr.LeastSquaresSolve(b);  // Minimize ||Ax-b||₂

// SVD (rank-deficient, singular matrices)
SVDecompositionSolver svd(A);
Vector<Real> x = svd.Solve(b);
int rank = svd.Rank();                // Matrix rank
Matrix<Real> null = svd.Nullspace();  // Null space basis
Real cond = 1.0/svd.inv_condition();  // Condition number

// Cholesky (symmetric positive definite - SPD)
CholeskySolver<Real> chol(A);         // A must be SPD!
Vector<Real> x = chol.Solve(b);
Real logdet = chol.logdet();          // log(det(A))

// Gauss-Jordan (small systems, need A⁻¹)
Matrix<Real> A_copy(A); Vector<Real> b_copy(b);
GaussJordanSolver<Real>::SolveInPlace(A_copy, b_copy);
// Result: A_copy = A⁻¹, b_copy = x
```

### Eigenvalues & Eigenvectors
```cpp
// Symmetric matrices (real eigenvalues)
SymmMatEigenSolverJacobi::Result result = SymmMatEigenSolverJacobi::Solve(A);
Vector<Real> eigenvals = result.eigenvalues;
Matrix<Real> eigenvecs = result.eigenvectors;
bool converged = result.converged;

// Access individual eigenvalues/eigenvectors
Real lambda_i = result.eigenvalues[i];
Vector<Real> v_i = result.eigenvectors.VectorFromColumn(i);

// General matrices: Power iteration method available in EigenSolverHelpers.h
```

---

## 📈 Calculus

### Derivatives
```cpp
RealFunction f = [](Real x) { return x*x; };

// First derivative f'(x₀) - different accuracy orders
Real df = Derivation::NDer1(f, x0);           // 1st order forward difference
Real df = Derivation::NDer2(f, x0);           // 2nd order central difference (recommended)
Real df = Derivation::NDer4(f, x0);           // 4th order central difference (high accuracy)

// Second derivative f''(x₀)
Real d2f = Derivation::NSecDer2(f, x0);       // 2nd order accuracy, 3 function evals
Real d2f = Derivation::NSecDer4(f, x0);       // 4th order accuracy, 5 function evals

// Third derivative f'''(x₀)
Real d3f = Derivation::NThirdDer2(f, x0);     // 2nd order accuracy, 4 function evals
Real d3f = Derivation::NThirdDer4(f, x0);     // 4th order accuracy, 6 function evals

// With error estimates
Real error;
Real df = Derivation::NDer2(f, x0, &error);

// Partial derivatives (multivariate) - see DerivationScalarFunction.h
```

### Integration (1D)
```cpp
RealFunction f = [](Real x) { return x*x; };

// Adaptive quadrature (recommended) - returns IntegrationResult
IntegrationResult result = IntegrateRomberg(f, a, b);     // Romberg (high accuracy)
IntegrationResult result = IntegrateSimpson(f, a, b);     // Simpson adaptive
IntegrationResult result = IntegrateTrap(f, a, b);        // Trapezoidal

// Access result
Real I = result.value;
bool converged = result.converged;
Real error = result.error;

// Gaussian quadrature - see GaussianQuadrature.h
Real I = GaussLegendre::Integrate(f, a, b, n);   // Legendre polynomials
Real I = GaussLaguerre::Integrate(f, 0, inf, n); // e⁻ˣ weight, [0,∞)
Real I = GaussHermite::Integrate(f, -inf, inf, n); // e⁻ˣ² weight, (-∞,∞)
```

### Integration (Multi-dimensional)
```cpp
// 2D integration
ScalarFunction<2> f = [](const VectorN<Real,2>& x) { ... };
IntegrationResult result = Integrate2D(f, method, ax, bx, ay, by, nx, ny);

// 3D integration
ScalarFunction<3> f = [](const VectorN<Real,3>& x) { ... };
IntegrationResult result = Integrate3D(f, method, ax,bx, ay,by, az,bz, nx,ny,nz);
```

---

## 🎯 Root Finding

### Polynomial Roots
```cpp
Polynom<Real> p{{1, -2, 1}};          // p(x) = x² - 2x + 1 (coef[0] = constant)

// Create from interpolation points
std::vector<Real> x{x0, x1, x2}, y{y0, y1, y2};
Polynom<Real> p = Polynom<Real>::FromValues(x, y);  // Lagrange interpolation

// Evaluate
Real val = p.Evaluate(x);
Real val = p(x);  // Operator overload

// Arithmetic
Polynom<Real> sum = p1 + p2;
Polynom<Real> prod = p1 * p2;
Polynom<Real> deriv = p.Derivative();

// Root finding: Use RootFinding algorithms with polynomial as function
RealFunction f = [&p](Real x) { return p(x); };
Real root = RootFinding::FindRootBrent(f, a, b, tol);
```

### Function Roots (1D)
```cpp
RealFunction f = [](Real x) { return x*x - 2; };
RealFunction df = [](Real x) { return 2*x; };    // Optional derivative

// Bracketing methods (no derivative needed)
Real root = RootFinding::FindRootBisection(f, a, b, tol);
Real root = RootFinding::FindRootBrent(f, a, b, tol);  // Recommended!

// Derivative-based (faster convergence)  
Real root = RootFinding::FindRootNewton(f, x0, x1, tol);  // Newton-Raphson
Real root = RootFinding::FindRootSecant(f, x0, x1, tol);  // Secant method
```

### Systems of Equations (N equations, N unknowns)
```cpp
// Newton-Raphson for F(x) = 0
VectorFunction<N> F = [...];
MatrixNM<Real,N,N> Jacobian(const VectorN<Real,N>& x) { ... }

VectorN<Real,N> x0 = {...};           // Initial guess
VectorN<Real,N> root = NewtonRaphson(F, Jacobian, x0);
```

---

## 🌊 Differential Equations

### Initial Value Problems (ODEs)
```cpp
// Define your ODE system (implements IODESystem interface)
class MyODESystem : public IODESystem {
public:
    int getDim() const override { return 2; }
    void derivs(Real t, const Vector<Real>& y, Vector<Real>& dydt) const override {
        dydt[0] = y[1];
        dydt[1] = -y[0];  // y'' + y = 0
    }
};

MyODESystem sys;
Vector<Real> y0{1.0, 0.0};  // Initial condition

// Fixed-step solvers
IODESystemStepCalculator* stepper = ...;  // Choose: EulerStep, RK4Step, etc.
ODESystemFixedStepSolver solver(sys, *stepper);
ODESystemSolution sol = solver.integrate(y0, t0, tend, numSteps);

// Adaptive solvers (recommended)
ODEAdaptiveIntegrator<DormandPrince5_Stepper> solver(sys);  // 5th order RK
ODESystemSolution sol = solver.integrate(y0, t0, tend, outputInterval, eps);

// Other adaptive steppers: CashKarp_Stepper, DormandPrince8_Stepper,
//                          BulirschStoer_Stepper, BulirschStoerRational_Stepper
```

---

## 📉 Optimization

### 1D Minimization
```cpp
RealFunction f = [](Real x) { return x*x + 3*x + 2; };

// Bracketing methods (no derivative)
Real xmin = Optimization::GoldenSection(f, a, b, tol);
Real xmin = Optimization::Brent(f, a, b, tol);     // Recommended!
```

### Multi-dimensional Minimization
```cpp
// Nelder-Mead simplex method (derivative-free)
class MyFunc : public IScalarFunction<N> {
public:
    Real operator()(const VectorN<Real,N>& x) const override { ... }
};

MyFunc func;
VectorN<Real,N> x0 = {...};           // Initial guess
Real delta = 0.1;                     // Initial simplex size

NelderMead<N> optimizer(maxIter, tol);
MultidimMinimizationResult result = optimizer.Minimize(func, x0, delta);

// Or use convenience function
MultidimMinimizationResult result = NelderMeadMinimize(func, x0, delta);

// Heuristic methods
SimulatedAnnealing<N> sa(config);
HeuristicOptimizationResult result = sa.Minimize(func, x0);

GeneticAlgorithm<N> ga(config);
HeuristicOptimizationResult result = ga.Minimize(func, x0);
```

---

## 🔊 Fourier & Signal Processing

### Fast Fourier Transform (FFT)
```cpp
Vector<Complex> signal(n);            // Input signal (n = power of 2)

// In-place transform
FFT::Transform(signal, 1);            // Forward FFT (isign=1)
FFT::Transform(signal, -1);           // Inverse FFT (isign=-1, need to normalize)

// Out-of-place
Vector<Complex> spectrum = FFT::TransformCopy(signal, 1);
Vector<Complex> restored = FFT::TransformCopy(spectrum, -1);
for (auto& val : restored) val /= n;  // Normalize inverse

// Real-valued FFT
Vector<Real> real_signal(n);
FFT::RealFFT(real_signal);            // Forward real FFT
FFT::RealIFFT(real_signal);           // Inverse real FFT

// Power spectrum
Vector<Real> power = FFT::PowerSpectrum(signal);
```

---

## 📊 Statistics

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

### Distributions
```cpp
// Normal distribution
Real pdf = NormalDistribution::PDF(x, mu, sigma);     // Probability density
Real cdf = NormalDistribution::CDF(x, mu, sigma);     // Cumulative distribution
Real x = NormalDistribution::InverseCDF(p, mu, sigma);// Quantile function

// Other distributions
Real p = BinomialDistribution::PMF(k, n, p);
Real p = PoissonDistribution::PMF(k, lambda);
Real p = ExponentialDistribution::PDF(x, lambda);
```

### Random Number Generation
```cpp
Ran ran(seed);                        // Initialize RNG with seed

// Uniform distributions
Real u = ran.doub();                  // Uniform [0,1)
int i = ran.int32();                  // Uniform 32-bit integer
int i = ran.int32(n);                 // Uniform [0, n)

// Continuous distributions
NormalDeviate norm(ran, mu, sigma);   // Normal(μ, σ)
Real x = norm.dev();

ExponentialDeviate exp(ran, lambda);  // Exponential(λ)
Real x = exp.dev();

// Discrete distributions
PoissonDeviate poisson(ran, lambda);  // Poisson(λ)
int k = poisson.dev();

BinomialDeviate binom(ran, n, p);     // Binomial(n, p)
int k = binom.dev();
```

---

## 🧮 Interpolation & Approximation

### Polynomial Interpolation
```cpp
std::vector<Real> x{x0, x1, x2, ...};
std::vector<Real> y{y0, y1, y2, ...};

// Lagrange interpolation via Newton's divided differences
Polynom<Real> p = Polynom<Real>::FromValues(x, y);
Real yi = p.Evaluate(xi);

// Or use InterpolatedFunction wrapper (see base/InterpolatedFunction.h)
```

### Least Squares Fitting
```cpp
// General least squares: minimize ||Ax - b||₂
QRSolver<Real> qr(A);
Vector<Real> x_fit = qr.LeastSquaresSolve(b);

// For polynomial fitting, build Vandermonde matrix manually:
// A[i][j] = x[i]^j for j=0..degree
```

---

## 📐 Geometry (2D/3D)

### 2D Primitives
```cpp
// Points
Point2Cartesian<Real> p(x, y);
Real dist = p.Distance(q);

// Lines
Line2D<Real> line(p1, p2);            // Line through two points
Real d = line.Distance(p);            // Point-to-line distance
bool intersects = line.Intersects(other_line);

// Polygons
Polygon2D<Real> poly{{p1, p2, p3, ...}};
Real area = poly.Area();
bool inside = poly.Contains(p);       // Point-in-polygon test
```

### 3D Primitives
```cpp
// Points
Point3Cartesian<Real> p(x, y, z);
Real dist = p.Distance(q);

// Lines
Line3D<Real> line(p1, p2);
Real d = line.Distance(p);

// Planes
Plane3D<Real> plane(p0, normal);
Real d = plane.Distance(p);
bool intersects = plane.Intersects(line);

// Triangles
Triangle3D<Real> tri(p1, p2, p3);
Real area = tri.Area();
Vector3Cartesian normal = tri.Normal();
```

---

## 🎲 Numerical Utilities

### Precision & Tolerances
```cpp
// Check if values are approximately equal
bool equal = Precision::IsEqual(a, b);            // Default: 1e-10
bool equal = Precision::IsEqual(a, b, tol);       // Custom tolerance

// Check if value is zero
bool is_zero = Precision::IsZero(x);              // |x| < ε
```

### Constants
```cpp
Real pi = Constants::PI;              // π
Real e = Constants::E;                // e (Euler's number)
Real phi = Constants::PHI;            // φ (golden ratio)
Real eps = Constants::EPSILON;        // Machine epsilon
```

---

## 💡 Quick Tips

**Linear Systems:** Use `CholeskySolver` for SPD matrices (fastest). Use `QRSolver` for ill-conditioned systems. Use `SVDecompositionSolver` for rank-deficient matrices.

**Integration:** Use `IntegrateRomberg` for smooth functions. Use Gaussian quadrature (GaussLegendre, etc.) for specific weights. For multi-dimensional, use adaptive methods or Monte Carlo.

**Root Finding:** Use `FindRootBrent` for 1D (no derivative). Use `FindRootNewton` when derivative available. Always bracket roots first.

**ODEs:** Use `ODEAdaptiveIntegrator` with `DormandPrince5_Stepper` for most problems. Use adaptive methods for stiff equations. Check stability for explicit methods.

**Optimization:** Use `Optimization::Brent` for 1D. Use `NelderMead` for multi-dimensional when derivatives unavailable.

**FFT:** Input size must be power of 2. Use `FFT::RealFFT` for real signals (2× faster). Zero-pad if needed.

---

## 📚 See Also

- **[README.md](../README.md)** — Installation, examples, getting started
- **[docs/](.)** — Complete documentation by module
- **[examples/](../examples/)** — Working code examples
- **[tests/](../tests/)** — Unit tests showing usage patterns

**Need more detail?** Check the full documentation or explore header files in `mml/`.

---

*MML v1.1 — Comprehensive C++ numerical computing library*
