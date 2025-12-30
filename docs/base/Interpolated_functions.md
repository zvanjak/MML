# Interpolated Functions

**Source:** `mml/base/InterpolatedFunction.h`

Interpolation creates continuous functions from discrete data points. MML provides several interpolation methods, each suited to different data characteristics and accuracy requirements.

## Overview

Available interpolation types:
- **Real functions** (IRealFunction): Linear, Polynomial, Rational, Barycentric, Cubic Spline
- **2D scalar functions**: Bilinear, Bicubic Spline
- **Parametric curves** (IParametricCurve<N>): Linear, Cubic Spline

## Quick Reference

| Method | Class | Data Type | Smoothness | Best For |
|--------|-------|-----------|------------|----------|
| Linear | `LinearInterpRealFunc` | 1D | C⁰ | Sparse data, noise tolerance |
| Polynomial | `PolynomInterpRealFunc` | 1D | C^(n-1) | Smooth data, low order |
| Rational | `RationalInterpRealFunc` | 1D | Varies | Near-pole behavior |
| Barycentric | `BarycentricRationalInterp` | 1D | C^∞ | Pole-free rational interp |
| Cubic Spline | `SplineInterpRealFunc` | 1D | C² | General purpose, smooth curves |
| Bilinear | `BilinearInterp2D` | 2D | C⁰ | Fast 2D, rectangular grids |
| Bicubic Spline | `BicubicSplineInterp2D` | 2D | C² | Smooth 2D surfaces |
| Linear Curve | `LinInterpParametricCurve<N>` | Parametric | C⁰ | Simple paths |
| Spline Curve | `SplineInterpParametricCurve<N>` | Parametric | C² | Smooth 3D paths |

---

## Real Function Interpolation

### LinearInterpRealFunc

Piecewise linear interpolation - connects adjacent data points with straight line segments.

**Constructor:**
```cpp
LinearInterpRealFunc(const Vector<Real>& x, const Vector<Real>& y, 
                     bool extrapolateOutsideOfRange = false);
```

**Parameters:**
- `x` - x-coordinates (must be monotonic - ascending or descending)
- `y` - y-coordinates
- `extrapolateOutsideOfRange` - if `false`, returns 0.0 outside [x_min, x_max]; if `true`, extrapolates linearly

**Properties:**
- Continuity: C⁰ (continuous, non-smooth at knots)
- Interpolation error: O(h²) where h is max spacing
- Complexity: O(log n) lookup + O(1) evaluation

**When to use:**
- ✅ Noisy or sparse data
- ✅ Data with discontinuous derivatives
- ✅ Fast evaluation required
- ❌ When smoothness is needed
- ❌ For highly oscillatory functions

**Example:**
```cpp
Vector<Real> x({ 0.0, 1.0, 2.0, 3.0, 4.0 });
Vector<Real> y({ 0.0, 1.0, 0.5, 2.0, 1.5 });

LinearInterpRealFunc f(x, y, false);  // no extrapolation

Real val = f(1.5);  // evaluates at x=1.5
```

### PolynomInterpRealFunc

Polynomial interpolation using Neville's algorithm with error estimation.

**Constructor:**
```cpp
PolynomInterpRealFunc(const Vector<Real>& x, const Vector<Real>& y, int m);
```

**Parameters:**
- `x` - x-coordinates (must be monotonic)
- `y` - y-coordinates  
- `m` - number of points used locally (polynomial order + 1)

**Properties:**
- Continuity: C^(m-1) (polynomial order m-1)
- Interpolation error: Available via `getLastErrorEst()`
- Complexity: O(log n) lookup + O(m²) evaluation
- Method: Neville's algorithm (numerically stable)

**When to use:**
- ✅ Smooth data with low polynomial order (m ≤ 4)
- ✅ When error estimates are needed
- ✅ Non-uniformly spaced data
- ❌ High-order polynomials (Runge's phenomenon)
- ❌ Data with noise

**Example:**
```cpp
Vector<Real> x({ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 });
Vector<Real> y({ 0.0, 0.8, 0.9, 0.1, -0.8, -1.0 });

PolynomInterpRealFunc f(x, y, 3);  // local cubic (order 3)

Real val = f(2.5);
Real err = f.getLastErrorEst();  // error estimate from last evaluation
```

**Warning - Runge's Phenomenon:**
High-order polynomial interpolation can produce wild oscillations between data points, especially near boundaries. Prefer m ≤ 4 or use splines instead.

### RationalInterpRealFunc

Diagonal rational function interpolation using the Bulirsch-Stoer algorithm. Useful for functions with poles or near-singular behavior.

**Constructor:**
```cpp
RationalInterpRealFunc(const Vector<Real>& x, const Vector<Real>& y, int m);
```

**Parameters:**
- `x` - x-coordinates (must be monotonic)
- `y` - y-coordinates
- `m` - number of points used locally for interpolation

**Properties:**
- Uses diagonal rational function through m nearest points
- Error estimation available via `getLastErrorEst()`
- Can interpolate functions with poles (unlike polynomials)
- Throws `RealFuncInterpRuntimeError` if evaluation hits a pole

**When to use:**
- ✅ Functions with poles or asymptotes
- ✅ Rational functions (ratios of polynomials)
- ✅ When polynomial interpolation oscillates excessively
- ❌ Data far from any pole-like behavior (use splines instead)
- ❌ When pole detection is not desired

**Example:**
```cpp
// Function with pole-like behavior: 1/(1 + x²)
Vector<Real> x({ -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0 });
Vector<Real> y({ 0.2, 0.5, 0.8, 1.0, 0.8, 0.5, 0.2 });

RationalInterpRealFunc f(x, y, 4);

Real val = f(0.25);
Real err = f.getLastErrorEst();  // error estimate

// Handle potential pole errors
try {
    Real val_near_pole = f(x_suspicious);
} catch (const RealFuncInterpRuntimeError& e) {
    // Pole detected in rational interpolation
}
```

**Reference:** Numerical Recipes §3.2 - Rational Function Interpolation

### BarycentricRationalInterp

Barycentric rational interpolation - a pole-free rational interpolation method that's numerically stable and guaranteed to have no poles within the data range.

**Constructor:**
```cpp
BarycentricRationalInterp(const Vector<Real>& x, const Vector<Real>& y, int d = 3);
```

**Parameters:**
- `x` - x-coordinates (must be monotonic)
- `y` - y-coordinates
- `d` - blending parameter (default 3), controls smoothness vs accuracy trade-off

**Properties:**
- **No poles guaranteed** within or near the data range
- Numerically stable (uses barycentric weights)
- Interpolates exactly at data points
- Smooth between data points
- O(n) evaluation complexity

**When to use:**
- ✅ When you need rational interpolation without pole risks
- ✅ Functions that vary smoothly
- ✅ When stability is paramount
- ✅ As an alternative to polynomial interpolation
- ❌ Functions that genuinely have poles (use `RationalInterpRealFunc`)

**Example:**
```cpp
Vector<Real> x({ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 });
Vector<Real> y({ 1.0, 1.5, 1.2, 2.0, 1.8, 2.5 });

BarycentricRationalInterp f(x, y, 3);  // d=3 (default)

Real val = f(2.5);  // guaranteed no pole

// Always exact at data points
for (int i = 0; i < x.size(); ++i) {
    assert(std::abs(f(x[i]) - y[i]) < 1e-12);
}
```

**Reference:** Berrut & Trefethen, "Barycentric Lagrange Interpolation" (SIAM Review, 2004)

### SplineInterpRealFunc

Cubic spline interpolation - piecewise cubic polynomials with continuous first and second derivatives. **The default choice for most interpolation tasks.**

**Constructor:**
```cpp
SplineInterpRealFunc(const Vector<Real>& x, const Vector<Real>& y,
                     Real yp1 = 1.e99, Real ypn = 1.e99);
```

**Parameters:**
- `x` - x-coordinates (must be monotonic)
- `y` - y-coordinates
- `yp1` - first derivative at x[0]; if ≥ 1e99, uses natural spline (y'' = 0)
- `ypn` - first derivative at x[n-1]; if ≥ 1e99, uses natural spline (y'' = 0)

**Properties:**
- Continuity: C² (continuous second derivative)
- Interpolation error: O(h⁴) for smooth functions
- Complexity: O(n) construction + O(log n) lookup + O(1) evaluation
- Storage: Stores second derivatives (_secDerY vector)

**Additional Methods:**
```cpp
Real Derivative(Real x) const;         // First derivative at x
Real SecondDerivative(Real x) const;   // Second derivative at x
Real Integrate(Real a, Real b) const;  // Exact integral from a to b
```

**When to use:**
- ✅ **DEFAULT CHOICE** for most interpolation tasks
- ✅ Smooth data requiring smooth interpolant
- ✅ Good balance of accuracy and stability
- ✅ Known or estimable endpoint derivatives
- ✅ When derivatives or integrals are needed
- ❌ Data with discontinuous derivatives
- ❌ Very large datasets (memory for second derivatives)

**Boundary Conditions:**
1. **Natural spline** (default): yp1 = ypn = 1e99 → second derivative = 0 at endpoints
2. **Clamped spline**: Specify yp1, ypn → first derivatives fixed at endpoints (more accurate if known)

**Example:**
```cpp
Vector<Real> x({ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 });
Vector<Real> y({ 0.0, 0.8, 0.9, 0.1, -0.8, -1.0, -0.5 });

// Natural spline (zero second derivative at boundaries)
SplineInterpRealFunc f_natural(x, y);

// Clamped spline (known derivatives at boundaries)
SplineInterpRealFunc f_clamped(x, y, 0.5, -0.3);

Real val = f_natural(2.5);

// Derivatives and integration
Real slope = f_natural.Derivative(2.5);
Real curvature = f_natural.SecondDerivative(2.5);
Real area = f_natural.Integrate(1.0, 4.0);
```

---

## 2D Interpolation

### BilinearInterp2D

Bilinear interpolation on a 2D rectangular grid. Fast and simple for smooth surfaces.

**Constructor:**
```cpp
BilinearInterp2D(const Vector<Real>& x, const Vector<Real>& y, 
                 const Matrix<Real>& z);
```

**Parameters:**
- `x` - x-coordinates of grid (must be monotonically increasing)
- `y` - y-coordinates of grid (must be monotonically increasing)
- `z` - function values at grid points; z(i,j) = f(x[i], y[j])

**Properties:**
- Continuity: C⁰ (continuous but not smooth at grid lines)
- Interpolation: Linear in each direction
- Complexity: O(log nx + log ny) lookup + O(1) evaluation
- Exact at grid points

**When to use:**
- ✅ Fast 2D interpolation
- ✅ Regular rectangular grids
- ✅ When smoothness is not critical
- ✅ Large grids (memory efficient)
- ❌ When C¹ or C² continuity is needed
- ❌ Highly curved surfaces

**Example:**
```cpp
// Create a 5x5 grid
Vector<Real> x({ 0.0, 1.0, 2.0, 3.0, 4.0 });
Vector<Real> y({ 0.0, 1.0, 2.0, 3.0, 4.0 });

// z values: z(i,j) = x[i] + y[j]
Matrix<Real> z(5, 5);
for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
        z(i, j) = x[i] + y[j];

BilinearInterp2D f(x, y, z);

// Evaluate at any point in the grid
Real val = f(1.5, 2.5);  // interpolated value

// Exact at grid points
Real exact = f(2.0, 3.0);  // = z(2,3) = 5.0
```

### BicubicSplineInterp2D

Bicubic spline interpolation on a 2D rectangular grid. Provides C² continuity for smooth surfaces.

**Constructor:**
```cpp
BicubicSplineInterp2D(const Vector<Real>& x, const Vector<Real>& y, 
                      const Matrix<Real>& z);
```

**Parameters:**
- `x` - x-coordinates of grid (must be monotonically increasing)
- `y` - y-coordinates of grid (must be monotonically increasing)
- `z` - function values at grid points; z(i,j) = f(x[i], y[j])

**Properties:**
- Continuity: C² (smooth surface with continuous first and second derivatives)
- Uses natural spline boundary conditions
- Precomputes spline coefficients for each row
- Complexity: O(log nx + log ny) lookup + O(ny) evaluation
- Exact at grid points

**When to use:**
- ✅ Smooth 2D surfaces
- ✅ When derivatives matter (gradients, curvature)
- ✅ Visualization of smooth terrain or fields
- ✅ Physical simulations requiring smooth potentials
- ❌ Very large grids (more memory than bilinear)
- ❌ When speed is critical (bilinear is faster)

**Example:**
```cpp
// Create grid for a smooth surface
Vector<Real> x({ 0.0, 1.0, 2.0, 3.0, 4.0 });
Vector<Real> y({ 0.0, 1.0, 2.0, 3.0, 4.0 });

// z values: smooth function z = sin(x) * cos(y)
Matrix<Real> z(5, 5);
for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
        z(i, j) = std::sin(x[i]) * std::cos(y[j]);

BicubicSplineInterp2D f(x, y, z);

// Evaluate - smooth between grid points
Real val = f(1.5, 2.5);

// Much smoother than bilinear for curved surfaces
```

---

## Parametric Curve Interpolation

Parametric curves are defined by position vectors **r**(t) = [x(t), y(t), z(t)]ᵀ in N-dimensional space.

### LinInterpParametricCurve<N>

Linear interpolation between N-dimensional points with arc-length parameterization.

**Constructor:**
```cpp
template<int N>
LinInterpParametricCurve(const Matrix<Real>& points, bool close = false);
```

**Parameters:**
- `points` - Matrix of curve points (numPoints × N); each row is one point
- `close` - if true, creates closed curve (wraps t parameter)

**Properties:**
- Parameterization: Arc-length normalized to [0,1]
- Continuity: C⁰ (continuous, piecewise linear)
- Parameter range: t ∈ [0,1]; wraps if closed

**Methods:**
```cpp
Real getMinT() const;        // returns 0.0
Real getMaxT() const;        // returns 1.0
VectorN<Real, N> operator()(Real t) const;  // evaluate at parameter t
```

**Example:**
```cpp
// 3D path through 5 points
Matrix<Real> points(5, 3, {
    0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 1.0, 0.0,
    0.0, 1.0, 0.0,
    0.0, 0.0, 1.0
});

LinInterpParametricCurve<3> curve(points, false);

VectorN<Real, 3> pos = curve(0.5);  // position at t=0.5
```

### SplineInterpParametricCurve<N>

Cubic spline interpolation for smooth N-dimensional curves with arc-length parameterization.

**Constructors:**
```cpp
template<int N>
SplineInterpParametricCurve(const Matrix<Real>& points, bool close = false);

template<int N>
SplineInterpParametricCurve(Real minT, Real maxT, const Matrix<Real>& points, 
                            bool close = false);
```

**Parameters:**
- `points` - Matrix of curve points (numPoints × N); each row is one point
- `close` - if true, creates closed curve; last point should NOT duplicate first
- `minT`, `maxT` - parameter range (default [0,1])

**Properties:**
- Each coordinate interpolated independently with cubic splines
- Continuity: C² in each coordinate
- Parameterization: Arc-length based, mapped to [minT, maxT]
- Endpoint derivatives: Estimated from neighboring points (≥4 points required)

**Methods:**
```cpp
Real getMinT() const;
Real getMaxT() const;
VectorN<Real, N> operator()(Real t) const;
```

**When to use:**
- ✅ Smooth 3D paths (camera paths, trajectories)
- ✅ Closed curves (loops)
- ✅ Animation paths
- ❌ Sharp corners (use linear interpolation)
- ❌ <4 points (needs neighbors for derivative estimation)

**Example:**
```cpp
// Smooth 3D path
Matrix<Real> points(10, 3, {
     0.0,  1.0,  5.0,
    -1.0,  2.0,  3.0, 
    -4.0,  3.0, -2.0, 
    20.0,  5.0, -2.0, 
     2.0,  3.0,  7.0, 
    -6.0, -3.0, -2.0, 
     3.0,  0.0,  5.0, 
     5.0, -2.0,  1.0, 
     6.0, -1.0, -2.0, 
    -2.0,  3.0,  0.0
});

SplineInterpParametricCurve<3> curve(points, false);

// Evaluate at 100 points for visualization
for (int i = 0; i <= 100; ++i) {
    Real t = i / 100.0;
    VectorN<Real, 3> pos = curve(t);
    // plot pos...
}

// Closed curve example (wraps around)
SplineInterpParametricCurve<2> closed_curve(loop_points, true);
```

---

## Algorithm Details

### Base Class: RealFunctionInterpolated

All real function interpolators derive from this base class.

**Internal Methods:**
- `locate(Real x)` - Binary search to find interval containing x (O(log n))
- `calcInterpValue(int startInd, Real x)` - Pure virtual, implemented by derived classes

**Data Storage:**
- Stores **copies** of input x and y arrays
- Sorts arrays if needed for monotonicity

### Spline Algorithm

Cubic spline construction solves a tridiagonal system for second derivatives:

1. **Input**: n points (xᵢ, yᵢ), boundary conditions (yp1, ypn)
2. **Solve**: Tridiagonal system for second derivatives y''ᵢ at each point
3. **Store**: Second derivative array (_secDerY)
4. **Evaluate**: For x ∈ [xⱼ, xⱼ₊₁], compute cubic polynomial using y, y' continuous at knots

**Complexity:**
- Construction: O(n) - tridiagonal solve
- Evaluation: O(log n) - binary search + O(1) cubic evaluation

**Reference:** Numerical Recipes §3.3

### Barycentric Rational Algorithm

The barycentric formula provides stable rational interpolation:

1. **Compute weights**: wⱼ = (-1)^j × C(d, j-k) for appropriate indices
2. **Evaluate**: r(x) = Σ wⱼyⱼ/(x-xⱼ) / Σ wⱼ/(x-xⱼ)

**Properties:**
- No poles in convex hull of data points
- Exact at interpolation nodes
- Numerically stable

**Reference:** Berrut & Trefethen, SIAM Review 2004

---

## Method Selection Guide

### Decision Tree

```
Do you have 1D data (y = f(x))?
├─ Yes → 
│  ├─ Is data noisy or sparse? → LinearInterpRealFunc
│  ├─ Function has poles/asymptotes? → RationalInterpRealFunc
│  ├─ Need pole-free rational? → BarycentricRationalInterp
│  ├─ Need smoothness? → SplineInterpRealFunc (DEFAULT)
│  ├─ Need error estimates? → PolynomInterpRealFunc (order ≤ 4)
│  └─ Data varies smoothly? → SplineInterpRealFunc
│
├─ Do you have 2D gridded data?
│  ├─ Need speed? → BilinearInterp2D
│  └─ Need smoothness? → BicubicSplineInterp2D
│
└─ Have parametric curve data (points in N-D space)?
   ├─ Need smooth curve? → SplineInterpParametricCurve<N>
   └─ Simple piecewise path? → LinInterpParametricCurve<N>
```

### Comparison Table - 1D Methods

| Criterion | Linear | Polynomial | Rational | Barycentric | Spline |
|-----------|--------|-----------|----------|-------------|--------|
| **Smoothness** | C⁰ | C^(m-1) | Varies | C^∞ | C² |
| **Accuracy** | O(h²) | O(h^m) | Variable | Good | O(h⁴) |
| **Stability** | Excellent | Poor (m>4) | Good | Excellent | Excellent |
| **Speed (eval)** | Fastest | Medium | Medium | O(n) | Fast |
| **Handles poles** | No | No | Yes | No (pole-free) | No |
| **Error estimate** | No | Yes | Yes | No | No |

### Comparison Table - 2D Methods

| Criterion | Bilinear | Bicubic Spline |
|-----------|----------|----------------|
| **Smoothness** | C⁰ | C² |
| **Accuracy** | O(h²) | O(h⁴) |
| **Speed** | Very fast | Fast |
| **Memory** | Low | Medium |
| **Best for** | Large grids, speed | Smooth surfaces |

### Common Pitfalls

**1. Runge's Phenomenon (Polynomial)**
```cpp
// ❌ DON'T: High-order polynomial on uniformly spaced data
PolynomInterpRealFunc bad(x, y, 12);  // oscillates wildly!

// ✅ DO: Use splines, rational, or low-order polynomial
SplineInterpRealFunc good(x, y);
BarycentricRationalInterp also_good(x, y);
PolynomInterpRealFunc ok(x, y, 3);
```

**2. Pole Detection (Rational)**
```cpp
// ✅ DO: Handle potential poles in rational interpolation
try {
    Real val = f_rational(x);
} catch (const RealFuncInterpRuntimeError& e) {
    // Handle pole case
}

// Or use BarycentricRationalInterp for guaranteed pole-free
BarycentricRationalInterp f_safe(x, y);
Real val = f_safe(x);  // never throws pole error
```

**3. Extrapolation**
```cpp
// ❌ DON'T: Extrapolate far outside data range
LinearInterpRealFunc f(x, y, true);  // extrapolation enabled
Real bad = f(x_max + 100);  // unreliable!

// ✅ DO: Only interpolate, or extrapolate cautiously
LinearInterpRealFunc f(x, y, false);  // returns 0 outside range
```

**4. Non-monotonic x values**
```cpp
// ❌ DON'T: Provide unsorted x values
Vector<Real> x({1.0, 3.0, 2.0, 4.0});  // not monotonic!

// ✅ DO: Ensure x is monotonically increasing or decreasing
Vector<Real> x({1.0, 2.0, 3.0, 4.0});
```

**5. 2D Grid Orientation**
```cpp
// ✅ DO: Ensure z(i,j) corresponds to f(x[i], y[j])
Matrix<Real> z(nx, ny);
for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
        z(i, j) = func(x[i], y[j]);  // correct indexing
```

---

## Performance Considerations

### Memory Usage

| Class | Storage | Notes |
|-------|---------|-------|
| LinearInterpRealFunc | 2n | x, y arrays |
| PolynomInterpRealFunc | 2n + 2m | x, y + workspace |
| RationalInterpRealFunc | 2n + 2m | x, y + workspace |
| BarycentricRationalInterp | 3n | x, y, weights |
| SplineInterpRealFunc | **3n** | x, y, y'' arrays |
| BilinearInterp2D | nx + ny + nx×ny | x, y, z |
| BicubicSplineInterp2D | nx + ny + 2×nx×ny | x, y, z, splines |
| Parametric curves | N×n + ... | points + splines |

### Evaluation Cost

```cpp
// For 1000 evaluations on n data points:

// Linear: ~1000 × log₂(n) comparisons + 1000 × 2 ops
LinearInterpRealFunc f(x, y);
// O(log n) each - fastest

// Barycentric: ~1000 × n ops (but very stable)
BarycentricRationalInterp f(x, y);
// O(n) each - slower but guaranteed pole-free

// Spline: ~1000 × log₂(n) + 1000 × 4 ops
SplineInterpRealFunc f(x, y);
// Fast evaluation - best speed/accuracy tradeoff

// 2D Bilinear: ~1000 × (log₂(nx) + log₂(ny)) + 1000 × 4 ops
BilinearInterp2D f(x, y, z);
// Very fast 2D lookup

// 2D Bicubic: ~1000 × (log₂(nx) + ny) ops
BicubicSplineInterp2D f(x, y, z);
// Slower than bilinear but much smoother
```

**Recommendation:** For repeated evaluations, splines offer best speed/accuracy tradeoff. For 2D, use bilinear for speed, bicubic for quality.

---

## Integration with Function Interface

All interpolation classes implement `IRealFunction`, `IScalarFunction2D`, or `IParametricCurve<N>`, so they work seamlessly with:

- **Derivation**: `NumericalDerivationAccurate`, `NumericalDerivationFast`
- **Integration**: `TrapezoidalIntegrator`, `SimpsonIntegrator`, `RombergIntegrator`
- **Visualization**: `Visualizer::VisualizeRealFunction`, `VisualizeParamCurve3D`
- **Serialization**: `Serializer::SaveRealFuncEquallySpaced`

**Example:**
```cpp
Vector<Real> x({0, 1, 2, 3, 4, 5});
Vector<Real> y({0, 1, 4, 9, 16, 25});

SplineInterpRealFunc f(x, y);

// Built-in derivative (exact for splines)
Real slope = f.Derivative(2.5);

// Built-in integral (exact for splines)
Real area = f.Integrate(0.0, 5.0);

// Or use numerical methods for other interpolation types
NumericalDerivationAccurate deriv;
Real slope2 = deriv(f, 2.5);

RombergIntegrator integrator;
Real area2 = integrator.integrate(f, 0.0, 5.0);

// Visualization
Visualizer::VisualizeRealFunction(f, "spline_plot", 0.0, 5.0, 100);
```

---

## Examples

### Example 1: Comparing 1D Interpolation Methods

```cpp
void Example_Compare_1D_Methods()
{
    Vector<Real> x({ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 });
    Vector<Real> y({ 0.0, 0.8, 0.9, 0.1, -0.8, -1.0 });

    LinearInterpRealFunc f_linear(x, y);
    PolynomInterpRealFunc f_poly(x, y, 3);
    RationalInterpRealFunc f_rational(x, y, 4);
    BarycentricRationalInterp f_bary(x, y);
    SplineInterpRealFunc f_spline(x, y);

    Real x_test = 2.5;
    std::cout << "Value at x = " << x_test << ":\n";
    std::cout << "  Linear:      " << f_linear(x_test) << "\n";
    std::cout << "  Polynomial:  " << f_poly(x_test) << "\n";
    std::cout << "  Rational:    " << f_rational(x_test) << "\n";
    std::cout << "  Barycentric: " << f_bary(x_test) << "\n";
    std::cout << "  Spline:      " << f_spline(x_test) << "\n";
}
```

### Example 2: 2D Surface Interpolation

```cpp
void Example_2D_Surface()
{
    // Create a 10x10 grid
    int n = 10;
    Vector<Real> x(n), y(n);
    for (int i = 0; i < n; ++i) {
        x[i] = i * 0.5;
        y[i] = i * 0.5;
    }

    // Sample a smooth function: z = sin(x) * cos(y)
    Matrix<Real> z(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            z(i, j) = std::sin(x[i]) * std::cos(y[j]);

    // Create interpolators
    BilinearInterp2D f_bilinear(x, y, z);
    BicubicSplineInterp2D f_bicubic(x, y, z);

    // Evaluate at intermediate point
    Real x_test = 1.25, y_test = 2.75;
    Real exact = std::sin(x_test) * std::cos(y_test);
    
    std::cout << "At (" << x_test << ", " << y_test << "):\n";
    std::cout << "  Exact:    " << exact << "\n";
    std::cout << "  Bilinear: " << f_bilinear(x_test, y_test) << "\n";
    std::cout << "  Bicubic:  " << f_bicubic(x_test, y_test) << "\n";
}
```

### Example 3: Spline Derivatives and Integration

```cpp
void Example_Spline_Calculus()
{
    // Data from sin(x) on [0, π]
    int n = 10;
    Vector<Real> x(n), y(n);
    for (int i = 0; i < n; ++i) {
        x[i] = i * Constants::PI / (n - 1);
        y[i] = std::sin(x[i]);
    }

    SplineInterpRealFunc f(x, y);

    // Derivative at π/4 (should be cos(π/4) ≈ 0.707)
    Real deriv = f.Derivative(Constants::PI / 4);
    std::cout << "f'(π/4) = " << deriv << " (exact: " << std::cos(Constants::PI/4) << ")\n";

    // Second derivative at π/2 (should be -sin(π/2) = -1)
    Real second = f.SecondDerivative(Constants::PI / 2);
    std::cout << "f''(π/2) = " << second << " (exact: -1)\n";

    // Integral from 0 to π (should be 2)
    Real integral = f.Integrate(0, Constants::PI);
    std::cout << "∫sin(x)dx = " << integral << " (exact: 2)\n";
}
```

### Example 4: Barycentric vs Rational for Stability

```cpp
void Example_Stability_Comparison()
{
    // Runge function: 1/(1 + 25x²) - notorious for polynomial oscillation
    int n = 11;
    Vector<Real> x(n), y(n);
    for (int i = 0; i < n; ++i) {
        x[i] = -1.0 + 2.0 * i / (n - 1);
        y[i] = 1.0 / (1.0 + 25.0 * x[i] * x[i]);
    }

    PolynomInterpRealFunc f_poly(x, y, n);  // Full polynomial - will oscillate!
    BarycentricRationalInterp f_bary(x, y); // Stable
    SplineInterpRealFunc f_spline(x, y);    // Also stable

    // Check at point where polynomial oscillates badly
    Real x_test = 0.9;
    Real exact = 1.0 / (1.0 + 25.0 * x_test * x_test);
    
    std::cout << "Runge function at x = " << x_test << ":\n";
    std::cout << "  Exact:       " << exact << "\n";
    std::cout << "  Polynomial:  " << f_poly(x_test) << " (oscillates!)\n";
    std::cout << "  Barycentric: " << f_bary(x_test) << " (stable)\n";
    std::cout << "  Spline:      " << f_spline(x_test) << " (stable)\n";
}
```

---

## Planned Features (TODO)

The following interpolation methods are **planned but not yet implemented**:

### Real Functions
- ❌ Hunt optimization for sequential lookups (O(1) vs O(log n))
- ❌ Monotonic spline (PCHIP) - preserves monotonicity of data
- ❌ Hermite spline - when derivative data is available

### N-dimensional
- ❌ Complete `InterpolatedSurface<N>` - currently a stub

---

## See Also

- [Functions.md](Functions.md) - Function interfaces (IRealFunction, IParametricCurve)
- [Derivation.md](Derivation.md) - Numerical derivatives of interpolated functions
- [Integration.md](Integration.md) - Numerical integration of interpolated functions
- [../base/README.md](../base/README.md) - Base data structures (Vector, Matrix)
- [../tools/Visualizers_README.md](../tools/Visualizers_README.md) - Visualization tools

---

## Runnable Examples

| Example | Description | Source |
|---------|-------------|--------|
| Linear Interpolation | LinearInterpRealFunc with equally spaced data | [`docs_demo_interpolated_functions.cpp`](../../src/docs_demos/docs_demo_interpolated_functions.cpp) |
| Spline Interpolation | SplineInterpRealFunc demonstration | [`docs_demo_interpolated_functions.cpp`](../../src/docs_demos/docs_demo_interpolated_functions.cpp) |
| Parametric Curves | SplineInterpParametricCurve<3> 3D path | [`docs_demo_interpolated_functions.cpp`](../../src/docs_demos/docs_demo_interpolated_functions.cpp) |
| 2D Bilinear | BilinearInterp2D surface interpolation | [`docs_demo_interpolated_functions.cpp`](../../src/docs_demos/docs_demo_interpolated_functions.cpp) |

**To run:** Build and execute `MML_DocsApp` target.

---

**References:**
- Press et al., *Numerical Recipes* (3rd ed.), Chapter 3: Interpolation and Extrapolation
- de Boor, *A Practical Guide to Splines*, Springer, 1978
- Berrut & Trefethen, "Barycentric Lagrange Interpolation", SIAM Review 46(3), 2004
