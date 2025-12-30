# Functions - Function Objects and Interfaces

**Files**: 
- `mml/interfaces/IFunction.h` - Function interfaces
- `mml/base/Function.h` - Function wrappers (lambda, function pointer, std::function)
- `mml/core/FunctionHelpers.h` - Function utilities and arithmetic operations

Function types and concrete wrappers for creating mathematical functions from lambdas, function pointers, and class members.

## Table of Contents
- [Overview](#overview)
- [Function Interfaces](#function-interfaces)
- [Creating Functions](#creating-functions)
- [Advanced Patterns](#advanced-patterns)
- [Examples](#examples)

---

## Overview

MML treats functions as first-class objects used across algorithms (derivation, integration, interpolation, ODEs, visualization). The library provides interfaces for different function types and convenient wrappers for practical use.

**Available Interfaces:**
- `IRealFunction` - f: ℝ → ℝ (1D real functions)
- `IScalarFunction<N>` - f: ℝⁿ → ℝ (scalar fields)
- `IVectorFunction<N>` - f: ℝⁿ → ℝⁿ (vector fields)
- `IVectorFunctionNM<N,M>` - f: ℝⁿ → ℝᵐ (general vector mappings)
- `IRealToVectorFunction<N>` - f: ℝ → ℝⁿ (parametric curves)
- `IParametricCurve<N>`, `IParametricSurface<N>` - geometric curves/surfaces
- `ITensorField2-5<N>` - tensor-valued fields (ranks 2-5)

**Key Features:**
- ✅ **Lambda support** - Create functions directly from lambda expressions
- ✅ **Function pointers** - Wrap existing C-style functions
- ✅ **std::function wrappers** - Support for complex capture scenarios
- ✅ **Class integration** - Use member functions and external class data
- ✅ **Automatic Jacobians** - Vector functions provide jacobian() method

---

## Function Interfaces

### IRealFunction (f: ℝ → ℝ)
```cpp
class IRealFunction : public IFunction<Real, Real> {
public:
    virtual Real operator()(Real x) const = 0;
};
```

**Use with**: 1D integration, derivation, root finding, 1D visualizers

### IScalarFunction\<N\> (f: ℝⁿ → ℝ)
```cpp
template<int N>
class IScalarFunction : public IFunction<Real, const VectorN<Real, N>&> {
public:
    virtual Real operator()(const VectorN<Real, N>& x) const = 0;
};
```

**Use with**: Multidimensional integration, gradient calculation, scalar field visualization

### IVectorFunction\<N\> (f: ℝⁿ → ℝⁿ)
```cpp
template<int N>
class IVectorFunction : public IFunction<VectorN<Real, N>, const VectorN<Real, N>&> {
public:
    virtual VectorN<Real, N> operator()(const VectorN<Real, N>& x) const = 0;
    virtual Real operator()(const VectorN<Real, N>& x, int component) const;
};
```

**Use with**: Jacobians, vector field operations, ODE systems, field visualization

### IVectorFunctionNM\<N,M\> (f: ℝⁿ → ℝᵐ)
```cpp
template<int N, int M>
class IVectorFunctionNM : public IFunction<VectorN<Real, M>, const VectorN<Real, N>&> {
public:
    virtual VectorN<Real, M> operator()(const VectorN<Real, N>& x) const = 0;
    virtual Real operator()(const VectorN<Real, N>& x, int component) const;
};
```

**Use with**: Coordinate transformations, general mappings, Jacobians (M×N matrices)

### IParametricCurve\<N\> (r: ℝ → ℝⁿ)
```cpp
template<int N>
class IParametricCurve : public IRealToVectorFunction<N> {
public:
    virtual VectorN<Real, N> operator()(Real t) const = 0;
    virtual Real getMinT() const = 0;
    virtual Real getMaxT() const = 0;
    
    std::vector<VectorN<Real, N>> GetTrace(double t1, double t2, int numPoints) const;
};
```

**Use with**: Curve visualization, arc length calculation, path integration

### IParametricSurface\<N\> (r: ℝ² → ℝⁿ)
```cpp
// Complex surface with variable w limits (depends on u)
template<int N>
class IParametricSurface : public IVectorFunctionNM<2, N> {
public:
    virtual VectorN<Real, N> operator()(Real u, Real w) const = 0;
    virtual VectorN<Real, N> operator()(const VectorN<Real, 2>& uv) const;
    
    virtual Real getMinU() const = 0;
    virtual Real getMaxU() const = 0;
    virtual Real getMinW(Real u) const = 0;  // w limits depend on u
    virtual Real getMaxW(Real u) const = 0;
};
```

### IParametricSurfaceRect\<N\> (r: ℝ² → ℝⁿ, rectangular domain)
```cpp
// Simple surface on rectangular domain
template<int N>
class IParametricSurfaceRect : public IVectorFunctionNM<2, N> {
public:
    virtual VectorN<Real, N> operator()(Real u, Real w) const = 0;
    virtual VectorN<Real, N> operator()(const VectorN<Real, 2>& coord) const;
    
    virtual Real getMinU() const = 0;
    virtual Real getMaxU() const = 0;
    virtual Real getMinW() const = 0;  // constant limits
    virtual Real getMaxW() const = 0;
};
```

**Use with**: Surface visualization, surface integrals, differential geometry

### ITensorField\<N\> (T: ℝⁿ → Tensor)
```cpp
template<int N>
class ITensorField2 : public IFunction<Tensor2<N>, const VectorN<Real, N>&> {
public:
    virtual Real Component(int i, int j, const VectorN<Real, N>& pos) const = 0;
    int getNumContravar() const;
    int getNumCovar() const;
};

// Similarly: ITensorField3, ITensorField4, ITensorField5
```

**Use with**: Stress/strain fields, metric tensors, curvature tensors

---

## Creating Functions

### Method 1: From Lambda Expressions (Simplest)
```cpp
// Real function
RealFunction f1([](Real x) { 
    return sin(x) * (1.0 + 0.5 * x * x); 
});

// Scalar field (3D)
ScalarFunction<3> potential([](const VectorN<Real,3>& x) {
    return 1.0 / x.NormL2();  // 1/r potential
});

// Vector field (2D)
VectorFunction<2> field([](const VectorN<Real,2>& x) {
    return VectorN<Real,2>{-x[1], x[0]};  // Rotation field
});

// Parametric curve (helix)
ParametricCurve<3> helix([](Real t) {
    return VectorN<Real,3>{cos(t), sin(t), t};
});

// Parametric surface (sphere on rectangular domain)
ParametricSurfaceRect<3> sphere([](Real u, Real v) {
    return VectorN<Real,3>{
        sin(u) * cos(v),
        sin(u) * sin(v),
        cos(u)
    };
}, 0, Constants::PI, 0, 2*Constants::PI);  // u: [0,π], v: [0,2π]

// Parametric surface with variable limits
ParametricSurface<3> cone([](Real u, Real v) {
    return VectorN<Real,3>{
        u * cos(v),
        u * sin(v),
        u
    };
}, 0, 1,                             // u: [0, 1]
   [](Real u) { return 0.0; },       // v_min(u) = 0
   [](Real u) { return 2*Constants::PI; });  // v_max(u) = 2π
```

### Method 2: From Function Pointers
```cpp
// Existing standalone function
Real my_function(Real x) {
    return exp(-x*x) * cos(x);
}

// Wrap it
RealFunction f2(my_function);

// Use with algorithms
Real integral = IntegrateSimpson(f2, 0.0, Constants::PI, 1e-6);
```

### Method 3: With std::function (Complex Captures)
```cpp
// Capture external state
double amplitude = 5.0;
double frequency = 2.0;

RealFunctionFromStdFunc wave(
    std::function<Real(Real)>([=](Real t) {
        return amplitude * sin(frequency * t);
    })
);

// Can modify captured values
amplitude = 10.0;  // But this won't affect 'wave' (captured by value)
```

### Method 4: From Class with Parameters
```cpp
// Custom function class with parameters
class TwoMassGravityPotential : public IScalarFunction<3> {
    Real m1, m2, G;
    VectorN<Real, 3> pos1, pos2;
    
public:
    TwoMassGravityPotential(Real m1_, Real m2_, Real G_,
                           const VectorN<Real,3>& p1,
                           const VectorN<Real,3>& p2)
        : m1(m1_), m2(m2_), G(G_), pos1(p1), pos2(p2) {}
    
    // Setters for parameter variation
    void SetM1(Real m) { m1 = m; }
    void SetM2(Real m) { m2 = m; }
    void SetPos1(const VectorN<Real,3>& p) { pos1 = p; }
    void SetPos2(const VectorN<Real,3>& p) { pos2 = p; }
    
    Real operator()(const VectorN<Real,3>& x) const override {
        Real r1 = (x - pos1).NormL2();
        Real r2 = (x - pos2).NormL2();
        return -G * (m1/r1 + m2/r2);
    }
};

// Usage
TwoMassGravityPotential potential(
    1000.0, 1000.0, 1.0,
    VectorN<Real,3>{10, 0, 0},
    VectorN<Real,3>{-10, 0, 0}
);

Real phi = potential(VectorN<Real,3>{0, 0, 0});

// Vary parameters
potential.SetM1(2000.0);
Real phi_new = potential(VectorN<Real,3>{0, 0, 0});
```

### Method 5: From External Class Member Functions
```cpp
// External class you can't modify
class PhysicsSimulation {
public:
    double temperature, pressure;
    
    double EnergyDensity(double x) const {
        return temperature * pressure * sin(x);
    }
};

// Wrap member function
PhysicsSimulation sim;
sim.temperature = 300.0;
sim.pressure = 101325.0;

RealFunctionFromStdFunc energy(
    std::function<Real(Real)>([&sim](Real x) {
        return sim.EnergyDensity(x);
    })
);

// Function updates when sim.temperature/pressure change
sim.temperature = 400.0;
Real e = energy(1.0);  // Uses new temperature
```

### Method 6: From Interpolated Data
```cpp
// Create function from sampled data
std::vector<Real> x_data = {0, 1, 2, 3, 4};
std::vector<Real> y_data = {0, 1, 4, 9, 16};

// Linear interpolation
LinearInterpRealFunc linear(x_data, y_data);
Real val = linear(2.5);  // Interpolates between points

// Spline interpolation
SplineInterpRealFunc spline(x_data, y_data);
Real smooth_val = spline(2.5);  // Smooth interpolation
```

---

## Advanced Patterns

### Jacobian Calculation
```cpp
// Vector function automatically computes Jacobian
VectorFunction<2> transform([](const VectorN<Real,2>& x) {
    return VectorN<Real,2>{
        x[0]*x[0] - x[1]*x[1],  // u = x² - y²
        2*x[0]*x[1]              // v = 2xy
    };
});

VectorN<Real,2> point{1.0, 2.0};

// Automatic Jacobian
MatrixNM<Real, 2, 2> J = transform.jacobian(point);
// J = [2x   -2y]   = [2   -4]
//     [2y    2x]     [4    2]
```

### Component Access
```cpp
VectorFunction<3> field([](const VectorN<Real,3>& x) {
    return VectorN<Real,3>{x[1], -x[0], x[2]};
});

VectorN<Real,3> pos{1, 2, 3};

// Full evaluation
VectorN<Real,3> value = field(pos);  // {2, -1, 3}

// Component-wise evaluation
Real comp_x = field(pos, 0);  // 2
Real comp_y = field(pos, 1);  // -1
Real comp_z = field(pos, 2);  // 3
```

### Parametric Curve Tracing
```cpp
// 3D helix
ParametricCurve<3> helix(0, 2*Constants::PI,
    [](Real t) {
        return VectorN<Real,3>{
            cos(t),
            sin(t),
            0.1 * t
        };
    }
);

// Generate trace points
auto trace = helix.GetTrace(0.0, 2*Constants::PI, 100);

// trace is std::vector<VectorN<Real,3>> with 100 points
```

---

## Function Helpers (FunctionHelpers.h)

### Taylor Series Expansion
```cpp
// Polynomial approximation around point 'a'
// f(x) ≈ f(a) + f'(a)(x-a) + f''(a)(x-a)²/2! + ...

RealFunction f([](Real x) { return sin(x); });

// 2nd order Taylor polynomial
PolynomRealFunc taylor2 = TaylorSeries2(f, 0.0);  // sin(x) ≈ x around 0

// 3rd order Taylor polynomial
PolynomRealFunc taylor3 = TaylorSeries3(f, 0.0);  // sin(x) ≈ x - x³/6
```

### Derivative Wrapper Classes
```cpp
// Create derivative as a function object
RealFunction f([](Real x) { return x*x*x; });  // f(x) = x³

// Different orders of accuracy (number = order of finite difference method)
RealFuncDerived1 df1(f);     // 1st order accuracy
RealFuncDerived2 df2(f);     // 2nd order accuracy
RealFuncDerived4 df4(f);     // 4th order accuracy
RealFuncDerived6 df6(f);     // 6th order accuracy
RealFuncDerived8 df8(f);     // 8th order accuracy

Real slope = df6(2.0);  // ≈ 12.0 (derivative of x³ at x=2)

// Custom step size
RealFuncDerived4 df4_fine(f, 1e-6);

// Second derivatives
RealFuncSecondDerived2 d2f2(f);  // 2nd order accuracy
RealFuncSecondDerived4 d2f4(f);  // 4th order accuracy
RealFuncSecondDerived6 d2f6(f);  // 6th order accuracy

Real curvature = d2f4(2.0);  // ≈ 12.0 (second derivative of x³ at x=2)
```

### Function Arithmetic Operations
```cpp
RealFunction f([](Real x) { return sin(x); });
RealFunction g([](Real x) { return cos(x); });

// Basic operations
RealFuncSum       h1(f, g);   // h(x) = sin(x) + cos(x)
RealFuncDiff      h2(f, g);   // h(x) = sin(x) - cos(x)
RealFuncProduct   h3(f, g);   // h(x) = sin(x) * cos(x)
RealFuncQuotient  h4(f, g);   // h(x) = sin(x) / cos(x) = tan(x)

// Composition: h(x) = f(g(x))
RealFuncCompose h5(f, g);     // h(x) = sin(cos(x))

// Scalar operations
RealFuncScale h6(f, 2.0);     // h(x) = 2*sin(x)
RealFuncShift h7(f, 1.0);     // h(x) = sin(x) + 1
RealFuncNegate h8(f);         // h(x) = -sin(x)

// Absolute and power
RealFuncAbs h9(f);            // h(x) = |sin(x)|
RealFuncPow h10(f, 2.0);      // h(x) = sin²(x)

// Comparison/distance
RealFuncAbsDiff h11(f, g);    // h(x) = |f(x) - g(x)|
RealFuncDiffSqr h12(f, g);    // h(x) = (f(x) - g(x))²  (for L² norms)
```

---

## Examples

### Example 1: Integration with Custom Function
```cpp
// Gaussian function
RealFunction gaussian([](Real x) {
    return exp(-x*x / 2.0) / sqrt(2.0 * Constants::PI);
});

// Integrate from -3σ to +3σ (≈ 0.9973)
Real probability = IntegrateSimpson(gaussian, -3.0, 3.0, 1e-8);

std::cout << "P(-3 < x < 3) = " << probability << std::endl;
```

### Example 2: Gradient of Scalar Field
```cpp
// 3D potential field: φ(x,y,z) = x² + y² + z²
ScalarFunction<3> potential([](const VectorN<Real,3>& x) {
    return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
});

VectorN<Real,3> point{1, 2, 3};

// Compute gradient: ∇φ = (2x, 2y, 2z)
auto grad = Gradient(potential, point);
// Result: {2, 4, 6}

std::cout << "∇φ at (1,2,3) = " << grad << std::endl;
```

### Example 3: Vector Field Divergence
```cpp
// 2D field: F(x,y) = (x, y) - radial field
VectorFunction<2> radial_field([](const VectorN<Real,2>& x) {
    return x;
});

VectorN<Real,2> point{3, 4};

// Divergence: ∇·F = ∂F_x/∂x + ∂F_y/∂y = 1 + 1 = 2
Real div = Divergence(radial_field, point);

std::cout << "∇·F = " << div << std::endl;  // 2
```

### Example 4: Parametric Surface with Physics
```cpp
// Temperature distribution on a sphere
class SphericalTemperature : public IParametricSurface<3> {
    Real radius;
    Real temp_equator, temp_pole;
    
public:
    SphericalTemperature(Real r, Real T_eq, Real T_pole)
        : radius(r), temp_equator(T_eq), temp_pole(T_pole) {}
    
    VectorN<Real,3> operator()(Real theta, Real phi) const override {
        // Position on sphere
        Real x = radius * sin(theta) * cos(phi);
        Real y = radius * sin(theta) * sin(phi);
        Real z = radius * cos(theta);
        
        // Temperature varies with latitude
        Real temp = temp_pole + (temp_equator - temp_pole) * sin(theta);
        
        // Color by temperature (could be used for visualization)
        return VectorN<Real,3>{x, y, z};
    }
    
    Real getMinX() const override { return 0; }  // theta: [0, π]
    Real getMaxX() const override { return Constants::PI; }
    Real getMinY() const override { return 0; }  // phi: [0, 2π]
    Real getMaxY() const override { return 2 * Constants::PI; }
};

// Earth-like temperature distribution
SphericalTemperature earth(6371.0, 300.0, 250.0);
```

### Example 5: ODE System from Vector Function
```cpp
// Lorenz attractor: dx/dt = σ(y-x), dy/dt = x(ρ-z)-y, dz/dt = xy-βz
Real sigma = 10.0, rho = 28.0, beta = 8.0/3.0;

VectorFunction<3> lorenz([=](const VectorN<Real,3>& state) {
    Real x = state[0], y = state[1], z = state[2];
    return VectorN<Real,3>{
        sigma * (y - x),
        x * (rho - z) - y,
        x * y - beta * z
    };
});

// Use with ODE solvers
VectorN<Real,3> initial{1, 1, 1};
auto solution = SolveODE_RK4(lorenz, initial, 0.0, 100.0, 0.01);
```

---

## Performance Notes

1. **Lambda vs Function Pointer**: Lambdas may inline better for simple functions
2. **std::function overhead**: Small overhead vs raw function pointers
3. **Jacobian calculation**: Uses finite differences (configurable step size)
4. **Component access**: Minimal overhead, evaluates full function once

---

## See Also
- [Interpolated_functions.md](../base/Interpolated_functions.md) - Create functions from data
- [Derivation.md](Derivation.md) - Numerical derivatives of functions
- [Integration.md](Integration.md) - 1D integration
- [Multidim_integration.md](Multidim_integration.md) - Multi-dimensional integration
- [Vector_field_operations.md](Vector_field_operations.md) - Gradient, divergence, curl
- [ODE_system.md](ODE_system.md) - Solving ODE systems
- [Curves_and_surfaces.md](Curves_and_surfaces.md) - Geometric applications

---

## Runnable Examples

| Example | Source File | Description |
|---------|------------|-------------|
| Functions Demo | [docs_demo_functions.cpp](../../src/docs_demos/docs_demo_functions.cpp) | All function creation patterns |

**Build and Run:**
```bash
cmake --build build --target MML_DocsApp
./build/src/docs_demos/Release/MML_DocsApp
```



