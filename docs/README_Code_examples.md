# 📝 MML Code Examples

Concise, copy-pasteable snippets showing the core MML API. Each section links to a full, compilable source file under [`../src/code_examples/`](../src/code_examples/).

> Back to the [main README](../README.md) • See also [Usage examples](README_Usage_examples.md) and the [Visualization suite](README_Visualization_suite.md).

---

## Contents

- [Vectors & Matrices](#vectors--matrices)
- [Linear Systems & Eigenvalues](#linear-systems--eigenvalues)
- [Polynomials & Algebra](#polynomials--algebra)
- [Defining Functions](#defining-functions)
- [Interpolation](#interpolation)
- [Numerical Derivatives](#numerical-derivatives)
- [Numerical Integration](#numerical-integration)
- [Root Finding](#root-finding)
- [Field Operations](#field-operations)
- [Coordinate Transformations](#coordinate-transformations)
- [Parametric Curves](#parametric-curves)
- [Path & Line Integrals](#path--line-integrals)
- [Function Analysis](#function-analysis)
- [Differential Equations](#differential-equations)
- [Dynamical Systems Analysis](#dynamical-systems-analysis)

---

## Vectors & Matrices

📄 *[View full source](../src/code_examples/readme02_vectors_matrices.cpp)*

```cpp
#include <MML.h>
using namespace MML;

// Real vectors and matrices
Vector<Real> vec1{1.5, -2.0, 0.5}, vec2{1.0, 1.0, -3.0};
Matrix<Real> mat_3x3{3, 3, {1.0,  2.0, -1.0,
                           -1.0,  5.0,  6.0,
                            3.0, -2.0,  1.0}};

// Vector and matrix arithmetic
Vector<Real> result = Real{2.0} * (vec1 + vec2) * mat_3x3 / vec1.NormL2();
std::cout << "Result: " << result << std::endl;

// Complex vectors and matrices
VectorComplex vec_cmplx{{Complex(1,1), Complex(-1,2)}};
MatrixComplex mat_cmplx{2, 2, {Complex(0.5,1),  Complex(-1,2),
                               Complex(-1,-2), Complex(-2,2)}};

// Matrix properties
std::cout << "IsOrthogonal:  " << MatrixAlg::IsOrthogonal(mat_3x3) << std::endl;
std::cout << "IsHermitian:   " << MatrixAlg::IsHermitian(mat_cmplx) << std::endl;

/* Expected OUTPUT:
    Result: [   -3.137858162,     3.922322703,    -8.629109946]
    IsOrthogonal:  0
    IsHermitian:   1
*/
```

---

## Linear Systems & Eigenvalues

📄 *[View full source](../src/code_examples/readme03_linear_systems.cpp)*

```cpp
// Define a linear system Ax = b
Matrix<Real> A{5, 5, {0.2,  4.1, -2.1, -7.4,  1.6,
                      1.6,  1.5, -1.1,  0.7,  5.0,
                     -3.8, -8.0,  9.6, -5.4, -7.8,
                      4.6, -8.2,  8.4,  0.4,  8.0,
                     -2.6,  2.9,  0.1, -9.6, -2.7}};
Vector<Real> b{1.1, 4.7, 0.1, 9.3, 0.4};

// Solve using LU decomposition
LUSolver<Real> luSolver(A);
Vector<Real> x = luSolver.Solve(b);

std::cout << "Solution: " << x << std::endl;
std::cout << "Verification A*x: " << (A * x) << std::endl;
std::cout << "Residual norm: " << (A*x - b).NormL2() << std::endl;

// Compute eigenvalues and eigenvectors
auto eigenResult = MatrixAlg::Eigensystem(A);

std::cout << "Eigenvalues:" << std::endl;
for (const auto& ev : eigenResult.eigenvalues)
    std::cout << "  λ = " << ev.real() << " + " << ev.imag() << "i" << std::endl;

// QR decomposition
QRSolver<Real> qr(A);
Matrix<Real> Q = qr.GetQ();
Matrix<Real> R = qr.GetR();
std::cout << "QR valid: " << Matrix<Real>::AreEqual(A, Q * R, 1e-10) << std::endl;
```

---

## Polynomials & Algebra

📄 *[View full source](../src/code_examples/readme09_polynomials.cpp)*

```cpp
// Create polynomials: p(x) = 2x³ - 3x² + x - 5
PolynomReal p{-5, 1, -3, 2};  // Coefficients: [constant, x, x², x³]

// Evaluate at a point
Real val = p(2.0);  // p(2) = 2(8) - 3(4) + 2 - 5 = 16 - 12 + 2 - 5 = 1
std::cout << "p(2) = " << val << std::endl;

// Polynomial arithmetic
PolynomReal q{1, 2};           // q(x) = 2x + 1
PolynomReal sum = p + q;       // Addition
PolynomReal prod = p * q;      // Multiplication
std::cout << "Degree of p*q: " << prod.degree() << std::endl;

// Calculus on polynomials
PolynomReal dp = p.derivative();        // p'(x) = 6x² - 6x + 1
PolynomReal ip = p.integral();          // ∫p(x)dx (constant = 0)
std::cout << "p'(x) at x=1: " << dp(1.0) << std::endl;

// Solve quadratic: x² - 5x + 6 = 0 → x = 2, 3
Complex r1, r2;
int numReal = SolveQuadratic(1.0, -5.0, 6.0, r1, r2);
std::cout << "Roots: " << r1.real() << ", " << r2.real() << std::endl;

// Solve cubic: x³ - 6x² + 11x - 6 = 0 → x = 1, 2, 3
Complex c1, c2, c3;
SolveCubic(1.0, -6.0, 11.0, -6.0, c1, c2, c3);
std::cout << "Cubic roots: " << c1.real() << ", " << c2.real()
          << ", " << c3.real() << std::endl;
```

---

## Defining Functions

📄 *[View full source](../src/code_examples/readme04_functions_interpolation.cpp)*

MML provides multiple ways to create function objects that can be derived, integrated, and analyzed:

```cpp
// =====================================================================
// CASE 1: From standalone function
// =====================================================================
double MyFunction(double x) { return sin(x) * (1.0 + 0.5*x*x); }
RealFunction f1(MyFunction);

// =====================================================================
// CASE 2: Direct lambda creation (most common and recommended)
// =====================================================================
RealFunction f2{[](Real x) { return sin(x) * (1.0 + 0.5*x*x); }};

// Different function types — all created with lambdas
ScalarFunction<3> scalar([](const VectorN<Real,3>& x) {
    return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];   // scalar field: r²
});

VectorFunction<3> vector([](const VectorN<Real,3>& x) -> VectorN<Real,3> {
    return {x[1]*x[2], x[0]*x[2], x[0]*x[1]};   // vector field
});

ParametricCurve<3> helix([](Real t) -> VectorN<Real,3> {
    return {cos(t), sin(t), 0.2*t};             // 3D parametric curve
});

ParametricSurface<3> torus([](Real u, Real v) -> VectorN<Real,3> {
    Real R = 3.0, r = 1.0;                      // parametric surface
    return {(R + r*cos(v))*cos(u), (R + r*cos(v))*sin(u), r*sin(v)};
});

// =====================================================================
// CASE 3a: From class with operator() — useful for stateful functions
// =====================================================================
class FunctionFromClassOperator {
    double _amplitude;
public:
    FunctionFromClassOperator(double amp) : _amplitude(amp) {}
    double operator()(double x) const { return _amplitude * std::sin(x); }
};

FunctionFromClassOperator obj(2.5);  // amplitude = 2.5
RealFunctionFromStdFunc f3(std::function<double(double)>{obj});

// =====================================================================
// CASE 3b: Class inheriting IRealFunction directly
// =====================================================================
class MyDerivedFunction : public IRealFunction {
    double _frequency;
public:
    MyDerivedFunction(double freq) : _frequency(freq) {}
    double operator()(double x) const override { return std::cos(_frequency * x); }
};
MyDerivedFunction f4(2.0);  // frequency = 2.0 → cos(2x)

// =====================================================================
// CASE 4: Wrapper for external/legacy class you can't modify
// =====================================================================
class ExternalComplexClass {
public:
    double ComputeValue(double x) const { return std::exp(-x*x); }
};

class ExternalClassWrapper : public IRealFunction {
    const ExternalComplexClass& _ref;
public:
    ExternalClassWrapper(const ExternalComplexClass& obj) : _ref(obj) {}
    double operator()(double x) const override { return _ref.ComputeValue(x); }
};

ExternalComplexClass externalObj;
ExternalClassWrapper f5(externalObj);  // Now usable with all MML algorithms!

// All these functions can now be derived, integrated, analyzed...
std::cout << "f2(1) = " << f2(1.0) << std::endl;           // 1.26147
std::cout << "f2'(1) = " << Derivation::NDer4(f2, 1.0) << std::endl;  // derivative
std::cout << "∫f2 = " << IntegrateSimpson(f2, 0, 1).value << std::endl;  // integral
```

---

## Interpolation

📄 *[View full source](../src/code_examples/readme04_functions_interpolation.cpp)*

Create smooth functions from discrete data points:

```cpp
// Data points (could be from experiment, simulation, file...)
Vector<Real> x_data{0, 2.5, 5.0, 7.5, 10.0};
Vector<Real> y_data{0, 1.5, 4.0, 9.5, 16.0};

// Interpolation methods create callable functions with value ownership
LinearInterpRealFunc   linear_interp(x_data, y_data);      // Piecewise linear
SplineInterpRealFunc   spline_interp(x_data, y_data);      // Cubic spline (smooth)
AkimaInterpRealFunc    akima_interp(x_data, y_data);       // Reduced noisy-data overshoot

// Evaluate anywhere in the interpolation range
Real x = 3.7;
std::cout << "Linear:     " << linear_interp(x) << std::endl;  // 2.70
std::cout << "Spline:     " << spline_interp(x) << std::endl;  // 2.41 (smoother)
std::cout << "Akima:      " << akima_interp(x) << std::endl;

// Detailed calls use one status/config contract across the 1D family
InterpolationConfig config;
config.extrapolation_policy = ExtrapolationPolicy::Clamp;
InterpolationResult detailed = akima_interp.EvaluateDetailed(12.0, config);
std::cout << detailed.value << " via " << detailed.algorithm_name << std::endl;
```

---

## Numerical Derivatives

📄 *[View full source](../src/code_examples/readme05_numerical_calculus.cpp)*

```cpp
// Compare derivative orders on f(x) = sin(x) at x = 1.0
// Analytical derivative: f'(1) = cos(1) ≈ 0.5403023058681398

RealFunction f{[](Real x) { return std::sin(x); }};
Real x = 1.0;
Real analytical = std::cos(1.0);

Real der1 = Derivation::NDer1(f, x);   // O(h) - forward difference
Real der2 = Derivation::NDer2(f, x);   // O(h²) - central difference
Real der4 = Derivation::NDer4(f, x);   // O(h⁴) - 5-point stencil
Real der6 = Derivation::NDer6(f, x);   // O(h⁶) - 7-point stencil
Real der8 = Derivation::NDer8(f, x);   // O(h⁸) - 9-point stencil

std::cout << std::scientific << std::setprecision(10);
std::cout << "NDer1 error: " << std::abs(der1 - analytical) << std::endl;  // ~1e-8
std::cout << "NDer2 error: " << std::abs(der2 - analytical) << std::endl;  // ~1e-11
std::cout << "NDer8 error: " << std::abs(der8 - analytical) << std::endl;  // ~1e-14
```

---

## Numerical Integration

📄 *[View full source](../src/code_examples/readme05_numerical_calculus.cpp)*

```cpp
// Integrate f(x) = sin(x) from 0 to π — Analytical result: 2.0
RealFunction f{[](Real x) { return std::sin(x); }};
Real a = 0.0, b = Constants::PI;

auto trap = IntegrateTrap(f, a, b, 1e-8);
auto simp = IntegrateSimpson(f, a, b, 1e-8);
auto romb = IntegrateRomberg(f, a, b, 1e-10);

std::cout << "Trapezoid: " << trap.value << " (error: " << trap.error_estimate << ")\n";
std::cout << "Simpson:   " << simp.value << " (error: " << simp.error_estimate << ")\n";
std::cout << "Romberg:   " << romb.value << " (error: " << romb.error_estimate << ")\n";

// Monte Carlo for high-dimensional integrals
ScalarFunction<3> volume_func([](const VectorN<Real, 3>& v) {
    return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
});
MonteCarloIntegrator<3> mc_integrator;
VectorN<Real, 3> lower{0, 0, 0}, upper{1, 1, 1};
auto mc_result = mc_integrator.integrate(volume_func, lower, upper,
                                          MonteCarloConfig().samples(100000));
std::cout << "MC estimate: " << mc_result.value
          << " +/- " << mc_result.error_estimate << std::endl;
```

---

## Root Finding

📄 *[View full source](../src/code_examples/readme10_root_finding.cpp)*

```cpp
// Find root of f(x) = x³ - 2x - 5 (has root near x ≈ 2.0945)
RealFunction f{[](Real x) { return x*x*x - 2*x - 5; }};

// Compare different methods
Real root_bisect = RootFinding::FindRootBisection(f, 2.0, 3.0, 1e-12);
Real root_brent  = RootFinding::FindRootBrent(f, 2.0, 3.0, 1e-12);
Real root_newton = RootFinding::FindRootNewton(f, 2.0, 3.0, 1e-12);
Real root_ridder = RootFinding::FindRootRidders(f, 2.0, 3.0, 1e-12);

std::cout << std::setprecision(15);
std::cout << "Bisection: " << root_bisect << std::endl;
std::cout << "Brent:     " << root_brent << std::endl;
std::cout << "Newton:    " << root_newton << std::endl;
std::cout << "Ridders:   " << root_ridder << std::endl;

// Get structured convergence diagnostics
RootFinding::RootFindingConfig config;
config.x_tolerance = 1e-14;
config.f_tolerance = 1e-14;
config.max_iterations = 100;

auto result = RootFinding::FindRootBrent(f, 2.0, 3.0, config);
std::cout << "Iterations: " << result.iterations_used << std::endl;
std::cout << "f(root) =   " << result.function_value << std::endl;

// Isolate and refine every detectable real root, including tangent roots
auto allRoots = RootFinding::FindAllRealRootsInInterval(
    [](Real x) { return (x - 1)*(x - 2)*(x - 2)*(x - 3); }, 0.0, 4.0);
std::cout << "Found " << allRoots.roots.size() << " unique roots" << std::endl;

// Solve F(x)=0 for a nonlinear system with a numerical Jacobian
auto system = [](const Vector<Real>& x) {
    return Vector<Real>({x[0]*x[0] + x[1]*x[1] - 1.0, x[0] - x[1]});
};
auto systemResult = RootFinding::SolveNonlinearSystemNewton(
    system, Vector<Real>({0.8, 0.4}));
```

> See [Root finding](algorithms/Root_finding.md) for scalar diagnostics, all-real-roots isolation, polynomial roots, and nonlinear system solvers.

---

## Field Operations

📄 *[View full source](../src/code_examples/readme06_field_operations.cpp)*

```cpp
// Scalar field: gravitational potential φ(x,y,z) = -1/r
ScalarFunction<3> potential([](const VectorN<Real, 3>& x) {
    return -1.0 / x.NormL2();
});
VectorN<Real, 3> pos{1.0, 2.0, 2.0};

// Gradient ∇φ (gives force direction)
auto grad = ScalarFieldOperations::GradientCart<3>(potential, pos);
std::cout << "Gradient at (1,2,2): " << grad << std::endl;

// Laplacian ∇²φ (zero outside mass for gravity!)
Real laplacian = ScalarFieldOperations::LaplacianCart<3>(potential, pos);
std::cout << "Laplacian at (1,2,2): " << laplacian << std::endl;

// Vector field: rotating velocity field v = (y, -x, z)
VectorFunction<3> velocity([](const VectorN<Real, 3>& x) -> VectorN<Real, 3> {
    return {x[1], -x[0], x[2]};
});

// Divergence ∇·v (compression/expansion rate)
Real div = VectorFieldOperations::DivCart<3>(velocity, pos);
std::cout << "Divergence at (1,2,2): " << div << std::endl;

// Curl ∇×v (rotation/vorticity)
auto curl = VectorFieldOperations::CurlCart(velocity, pos);
std::cout << "Curl at (1,2,2): " << curl << std::endl;
```

---

## Coordinate Transformations

📄 *[View full source](../src/code_examples/readme11_coord_transforms.cpp)*

```cpp
// Cartesian point
Vector3Cartesian cart_pos{1.0, 1.0, 1.0};

// Convert to Spherical (r, θ, φ) - Math/ISO convention
// θ = polar angle from z-axis, φ = azimuthal angle in xy-plane
Vector3Spherical sph_pos = CoordTransfCartToSpher.transf(cart_pos);
std::cout << "Spherical: r=" << sph_pos[0] << ", θ=" << sph_pos[1]
          << ", φ=" << sph_pos[2] << std::endl;

// Convert to Cylindrical (r, φ, z)
Vector3Cylindrical cyl_pos = CoordTransfCartToCyl.transf(cart_pos);
std::cout << "Cylindrical: r=" << cyl_pos[0] << ", φ=" << cyl_pos[1]
          << ", z=" << cyl_pos[2] << std::endl;

// Convert back to Cartesian
Vector3Cartesian back = CoordTransfSpherToCart.transf(sph_pos);
std::cout << "Back to Cartesian: " << back << std::endl;

// Field operations in spherical coordinates: φ = -1/r
ScalarFunction<3> pot_spher([](const VectorN<Real, 3>& x) { return -1.0/x[0]; });

Vector3Spherical test_sph{2.0, Constants::PI/4, Constants::PI/4};
auto grad_spher = ScalarFieldOperations::GradientSpher(pot_spher, test_sph);
std::cout << "  ∂φ/∂r = " << grad_spher[0] << " (analytical: 1/r² = 0.25)" << std::endl;

// Covariant vector transformation (transform gradient to Cartesian coords)
Vector3Cartesian test_cart = CoordTransfSpherToCart.transf(test_sph);
auto force_cart = CoordTransfSpherToCart.transfVecCovariant(grad_spher, test_cart);
std::cout << "Force vector in Cartesian: " << force_cart << std::endl;
```

---

## Parametric Curves

📄 *[View full source](../src/code_examples/readme12_parametric_curves.cpp)*

```cpp
// Define a 3D helix: r(t) = (cos(t), sin(t), 0.2t)
Curves::CurveCartesian3D helix([](Real t) -> VectorN<Real, 3> {
    return {cos(t), sin(t), 0.2*t};
});

Real t = Constants::PI / 4;

// Curve properties at parameter t
auto pos = helix(t);                    // Position on curve
auto tangent = helix.getTangent(t);     // Tangent vector dr/dt
auto unit_tan = helix.getTangentUnit(t);// Unit tangent T
auto normal = helix.getNormal(t);       // Normal vector (acceleration)
auto binormal = helix.getBinormal(t);   // Binormal B = T × N

std::cout << std::setprecision(6);
std::cout << "Position:       " << pos << std::endl;
std::cout << "Unit tangent:   " << unit_tan << std::endl;
std::cout << "Binormal:       " << binormal << std::endl;

// Curvature κ (Frenet-Serret apparatus)
Real curvature = helix.getCurvature(t);
std::cout << "Curvature κ:    " << curvature << std::endl;

// Predefined curves
Curves::LemniscateCurve lemniscate;       // Figure-eight curve
Curves::ToroidalSpiralCurve torus(5, 2);  // Spiral on torus
Curves::Circle3DXZCurve circle(3.0);      // Circle in XZ plane
```

---

## Path & Line Integrals

📄 *[View full source](../src/code_examples/readme13_path_integrals.cpp)*

```cpp
// Line integral of vector field along a curve: ∫ F·dr, F = (y, -x, 0)
VectorFunction<3> F([](const VectorN<Real, 3>& p) -> VectorN<Real, 3> {
    return {p[1], -p[0], 0.0};
});

// Unit circle in XY plane: r(t) = (cos(t), sin(t), 0), t ∈ [0, 2π]
ParametricCurve<3> circle([](Real t) -> VectorN<Real, 3> {
    return {cos(t), sin(t), 0.0};
});

// Work integral: W = ∫₀^{2π} F(r(t)) · r'(t) dt
Real work = PathIntegration::LineIntegral(F, circle, 0.0, 2*Constants::PI, 1e-8);
std::cout << "Work ∮ F·dr: " << work << " (analytical -2π = " << -2*Constants::PI << ")\n";

// Arc length computation: ∫ ds
Real arc_length = PathIntegration::ParametricCurveLength<3>(circle, 0.0, 2*Constants::PI);
std::cout << "Arc length: " << arc_length << " (analytical 2π = " << 2*Constants::PI << ")\n";
```

---

## Function Analysis

📄 *[View full source](../src/code_examples/readme14_function_analysis.cpp)*

```cpp
// Analyze function behavior over an interval
RealFunction f{[](Real x) { return x*x*x - 3*x + 1; }};

RealFunctionAnalyzer analyzer(f, "x³ - 3x + 1");
analyzer.PrintIntervalAnalysis(-3.0, 3.0, 100, 1e-6);
/*  Output:
    f(x) = x³ - 3x + 1 - Analysis in interval [-3.00, 3.00]:
      Defined    : yes
      Continuous : yes
      Monotonic  : no
      Min        : -0.999...
      Max        : 2.999...
*/

// Find roots using root bracket search + bisection
Vector<Real> xb1, xb2;
int numBrackets = RootFinding::FindRootBrackets(f, -3.0, 3.0, 100, xb1, xb2);
std::cout << std::setprecision(10);
for (int i = 0; i < numBrackets; i++) {
    Real root = RootFinding::FindRootBisection(f, xb1[i], xb2[i], 1e-10);
    std::cout << "  x = " << root << ", f(x) = " << f(root) << std::endl;
}

// Analyze a function with discontinuity
RealFunctionFromStdFunc step([](Real x) -> Real {
    if (x < 0) return 0.0;
    else if (x > 0) return 1.0;
    else return 0.5;
});
RealFunctionAnalyzer step_analyzer(step, "step(x)");
step_analyzer.PrintIntervalAnalysis(-2.0, 2.0, 100, 1e-6);  // Detects discontinuity at x = 0
```

---

## Differential Equations

📄 *[View full source](../src/code_examples/readme07_ode_solvers.cpp)*

```cpp
// Simple harmonic oscillator: d²x/dt² = -ω²x  →  dx/dt = v, dv/dt = -ω²x
ODESystem system(2, [](Real t, const Vector<Real>& y, Vector<Real>& dydt) {
    const Real omega = 2.0;
    dydt[0] = y[1];                    // dx/dt = v
    dydt[1] = -omega*omega * y[0];     // dv/dt = -ω²x
});

Vector<Real> initial_cond{1.0, 0.0};  // x(0) = 1, v(0) = 0

// Solve with RK4 stepper
RungeKutta4_StepCalculator stepper;
ODESystemFixedStepSolver solver(system, stepper);
ODESystemSolution solution = solver.integrate(initial_cond, 0.0, 10.0, 1000);

int last = solution.size() - 1;
std::cout << "Final position: " << solution.getXValue(last, 0) << std::endl;
std::cout << "Final velocity: " << solution.getXValue(last, 1) << std::endl;

// Analytical solution at t=10: x(t) = cos(ωt), v(t) = -ω sin(ωt)
const Real omega = 2.0, t_final = 10.0;
std::cout << "Analytical x:   " << std::cos(omega * t_final) << std::endl;
std::cout << "Analytical v:   " << -omega * std::sin(omega * t_final) << std::endl;
```

---

## Dynamical Systems Analysis

📄 *[View full source](../src/code_examples/readme15_dynamical_systems.cpp)*

```cpp
// THE LORENZ SYSTEM — the butterfly effect in action
// dx/dt = σ(y - x), dy/dt = x(ρ - z) - y, dz/dt = xy - βz
using namespace MML::Systems;

LorenzSystem lorenz(10.0, 28.0, 8.0/3.0);  // Classic chaotic parameters

std::cout << "Dissipative: " << (lorenz.isDissipative() ? "yes" : "no") << std::endl;
std::cout << "Flow divergence: " << lorenz.getDivergence() << std::endl;

// FIXED POINT ANALYSIS — find equilibria and classify stability
std::vector<Vector<Real>> guesses = {
    Vector<Real>{0.0, 0.0, 0.0}, Vector<Real>{8.0, 8.0, 27.0}, Vector<Real>{-8.0, -8.0, 27.0}
};
auto fixedPoints = FixedPointFinder::FindMultiple(lorenz, guesses);
for (const auto& fp : fixedPoints)
    std::cout << "  Fixed point type: " << ToString(fp.type)
              << ", stable: " << (fp.isStable ? "yes" : "NO") << std::endl;

// LYAPUNOV EXPONENTS — quantify chaos via sensitivity to initial conditions
Vector<Real> x0 = lorenz.getDefaultInitialCondition();
auto lyap = LyapunovAnalyzer::Compute(lorenz, x0, 500.0, 1.0, 0.01);
std::cout << "Max exponent: " << lyap.maxExponent
          << (lyap.maxExponent > 0 ? " (CHAOS!)" : "") << std::endl;
std::cout << "Kaplan-Yorke dimension: " << lyap.kaplanYorkeDimension << std::endl;

// BIFURCATION ANALYSIS — route to chaos via parameter sweeps
auto bifurcation = BifurcationAnalyzer::Sweep(
    lorenz, 1 /*rho*/, 20.0, 30.0, 6, Vector<Real>{1.0, 1.0, 1.0}, 2, 50.0, 20.0, 0.01);
// Reveals: rho < 24 periodic → period-doubling cascade → rho ≈ 24.74 chaos onset
```

> See [Dynamical Systems](systems/) for fixed-point classification, Lyapunov spectra, and bifurcation analysis.