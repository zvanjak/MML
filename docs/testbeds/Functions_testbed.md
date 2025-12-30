# Functions Test Bed

## Overview

The Functions Test Bed provides a comprehensive collection of mathematical functions for validating numerical algorithms including differentiation, integration, optimization, and vector calculus operations. Each test function includes analytically-derived reference values for derivatives, integrals, and other properties.

The test bed is organized into three main categories:

| Test Bed | Description | Count |
|----------|-------------|-------|
| **RealFunctionsTestBed** | f: ℝ → ℝ (single-variable functions) | 60+ functions |
| **ScalarFunctionsTestBed** | f: ℝⁿ → ℝ (multi-variable scalar functions) | 15 functions |
| **VectorFunctionsTestBed** | F: ℝⁿ → ℝⁿ (vector fields) | 5 functions |

## Architecture

The test bed is organized across header files in `test_data/`:

| File | Purpose |
|------|---------|
| `real_functions_test_bed.h` | Single-variable functions with derivatives and integrals |
| `scalar_functions_test_bed.h` | Multi-variable scalar functions with gradients |
| `vector_functions_test_bed.h` | Vector fields with Jacobians and curl/divergence |

---

# Real Functions Test Bed

## Data Structures

### TestFunctionReal (Complete Function Data)

The most comprehensive test container, including function through third derivative and integral:

```cpp
struct TestFunctionReal {
    std::string _funcName;
    
    std::shared_ptr<IInterval> _intervalDef;   // Mathematical domain
    std::shared_ptr<IInterval> _intervalTest;  // Testing interval (subset of domain)
    
    RealFunction _func;           // f(x)
    RealFunction _funcDerived;    // f'(x)
    RealFunction _funcSecDer;     // f''(x)
    RealFunction _funcThirdDer;   // f'''(x)
    RealFunction _funcIntegrated; // ∫f(x)dx
    
    std::string _funcExpr;        // Expression strings for display
    std::string _funcDerivedExpr;
    std::string _funcSecDerExpr;
    std::string _funcThirdDerExpr;
    std::string _funcIntegratedExpr;
};
```

### TestFunctionRealWithDerivation

Lighter structure for derivative testing only:

```cpp
struct TestFunctionRealWithDerivation {
    std::string _funcName;
    std::shared_ptr<IInterval> _intervalDef, _intervalTest;
    RealFunctionFromStdFunc _func;
    RealFunctionFromStdFunc _funcDerived;
    std::string _funcExpr, _funcDerivedExpr;
};
```

### TestFunctionRealWithIntegral

Lighter structure for integration testing only:

```cpp
struct TestFunctionRealWithIntegral {
    std::string _funcName;
    std::shared_ptr<IInterval> _intervalDef, _intervalTest;
    RealFunctionFromStdFunc _func;
    RealFunctionFromStdFunc _funcIntegrated;
    std::string _funcExpr, _funcIntegratedExpr;
};
```

## Standard Functions (19 functions)

Complete functions with derivatives through third order and antiderivatives:

| Name | Function | Domain | Test Interval |
|------|----------|--------|---------------|
| Sin | sin(x) | ℝ | [-2π, 2π] |
| Cos | cos(x) | ℝ | [-2π, 2π] |
| Tan | tan(x) | ℝ \ {π/2 + nπ} | [-5π, 5π] with holes |
| Sinh | sinh(x) | ℝ | [-10, 10] |
| Cosh | cosh(x) | ℝ | [-10, 10] |
| Tanh | tanh(x) | ℝ | [-10, 10] |
| Sqrt | √x | (0, ∞) | (0, 10⁶] |
| x² | x² | ℝ | [-100, 100] |
| x³ | x³ | ℝ | [-100, 100] |
| x⁴ | x⁴ | ℝ | [-100, 100] |
| x⁵ | x⁵ | ℝ | [-100, 100] |
| Ln | ln(x) | (0, ∞) | (0, 1000] |
| Exp | eˣ | ℝ | (-20, 20) |
| Asin | arcsin(x) | (-1, 1) | (-1, 1) |
| Acos | arccos(x) | (-1, 1) | (-1, 1) |
| Atan | arctan(x) | ℝ | [-100, 100] |
| Asinh | arcsinh(x) | ℝ | [-10, 10] |
| Acosh | arccosh(x) | (-∞,-1) ∪ (1,∞) | Intervals avoiding ±1 |
| Atanh | arctanh(x) | (-1, 1) ∪ exterior | Multiple intervals |

**Purpose:** Core mathematical functions for baseline testing of differentiation and integration algorithms.

## Derivative Test Functions (11 functions)

Functions specifically chosen for derivative testing:

| Name | f(x) | f'(x) | Test Interval |
|------|------|-------|---------------|
| TestDer1 | sin(x) | cos(x) | (-20, 20) |
| TestDer2 | x³ | 3x² | (-20, 20) |
| TestDer3 | tanh(x) | 1/cosh²(x) | (-10, 10) |
| Gaussian | e^(-x²) | -2xe^(-x²) | [-5, 5] |
| Logistic | 1/(1+e^(-x)) | e^(-x)/(1+e^(-x))² | [-10, 10] |
| ErrorFunction | erf(x) | (2/√π)e^(-x²) | [-3, 3] |
| ReciprocalQuadratic | 1/(1+x²) | -2x/(1+x²)² | [-10, 10] |
| PolyMixedPowers | x⁵-3x³+2x | 5x⁴-9x²+2 | [-10, 10] |
| SinCosComposite | sin(x)cos(x) | cos²(x)-sin²(x) | [-2π, 2π] |
| ExpPoly | x²eˣ | eˣ(x²+2x) | [-10, 10] |
| RationalFunc | (x²+1)/(x³+2) | [complex] | [-10, 10] |

**Purpose:** Tests numerical differentiation with various function types including transcendental, rational, and composite functions.

## Integration Test Functions (16 functions)

Functions with known antiderivatives:

| Name | f(x) | ∫f(x)dx | Test Interval |
|------|------|---------|---------------|
| TestInt1 | x²(x²-2)sin(x) | 4x(x²-7)sin(x)-(x⁴-14x²+28)cos(x) | (0, 5) |
| TestInt2 | 1/(4+x²) | ½atan(x/2) | (-20, 20) |
| XSquared | x² | x³/3 | (-10, 10) |
| Reciprocal | 1/x | ln|x| | (0.1, 10) |
| ExpDecay | e^(-2x) | -½e^(-2x) | [-10, 10] |
| Gaussian | e^(-x²) | ½√π·erf(x) | [-3, 3] |
| Logistic | 1/(1+e^(-x)) | x - ln(1+e^(-x)) | [-10, 10] |
| RationalXOverXSquaredPlus1 | x/(x²+1) | ½ln(x²+1) | [-10, 10] |
| Arctangent | atan(x) | x·atan(x) - ½ln(1+x²) | [-10, 10] |
| SineSquared | sin²(x) | x/2 - ¼sin(2x) | [-10, 10] |
| CosineCubed | cos³(x) | sin(x)/3 + sin(3x)/9 | [-10, 10] |
| ReciprocalQuadratic | 1/(x²+4) | ½atan(x/2) | [-10, 10] |
| ExpTimesX | xeˣ | (x-1)eˣ | [-10, 10] |
| LogComposite | ln(1+x²) | x·ln(1+x²) - 2x + 2·atan(x) | [-10, 10] |
| Tanh | tanh(x) | ln(cosh(x)) | [-10, 10] |
| ExpQuadratic | e^(x²) | ½√π·erf(x) | [-3, 3] |

**Purpose:** Validates numerical integration algorithms with various function types.

## Numerically Challenging Functions (8 extended)

Functions specifically chosen for their challenging numerical properties:

### 1. Runge Function
```
f(x) = 1/(1 + 25x²)
f'(x) = -50x/(1 + 25x²)²
```
**Properties:** Classic example of Runge's phenomenon - polynomial interpolation fails at edges. Peaked function that challenges interpolation and quadrature.

### 2. Steep Gradient (k=10)
```
f(x) = tanh(10x)
f'(x) = 10/cosh²(10x)
```
**Properties:** Rapid transition region tests adaptive methods and derivative approximations near steep gradients.

### 3. Peaked Gaussian (ε=0.01)
```
f(x) = exp(-x²/0.01)
f'(x) = -200x·exp(-x²/0.01)
```
**Properties:** Extremely narrow peak (width ~0.1) challenges quadrature rules to detect and resolve the feature.

### 4. Witch of Agnesi
```
f(x) = 1/(1 + x²)
f'(x) = -2x/(1 + x²)²
```
**Properties:** Classical benchmark function, well-behaved but with asymptotic tails.

### 5. Oscillatory Decay
```
f(x) = exp(-x)·sin(10x)
f'(x) = exp(-x)·(10·cos(10x) - sin(10x))
```
**Properties:** Combines exponential decay with high-frequency oscillation - challenges many algorithms.

### 6. Near-Singular Power
```
f(x) = x^0.1
f'(x) = 0.1·x^(-0.9)
```
**Properties:** Infinite derivative at x=0 tests boundary handling.

### 7. Sinc Function
```
f(x) = sin(x)/x  (with f(0) = 1)
f'(x) = (x·cos(x) - sin(x))/x²
```
**Properties:** Removable singularity at origin, oscillatory decay.

### 8. Smooth Bump
```
f(x) = exp(-1/(1-x²))  for |x| < 1
```
**Properties:** C^∞ function with compact support - infinitely differentiable but very flat near boundaries.

## Challenging Integration Functions (6 extended)

### 1. Fresnel Sine
```
f(x) = sin(x²)
∫f(x)dx = √(π/2)·S(√(2/π)·x)  [Fresnel S function]
```
**Properties:** Highly oscillatory with increasing frequency - standard quadrature struggles.

### 2. Fresnel Cosine
```
f(x) = cos(x²)
∫f(x)dx = √(π/2)·C(√(2/π)·x)  [Fresnel C function]
```
**Properties:** Same challenges as Fresnel sine.

### 3. Runge Integral
```
f(x) = 1/(1+25x²)
∫f(x)dx = (1/5)·atan(5x)
```
**Properties:** Integration of peaked function.

### 4. Peaked Gaussian Integral
```
f(x) = exp(-x²/0.01)
∫f(x)dx = 0.5·√(0.01π)·erf(10x)
```
**Properties:** Narrow peak requires adaptive refinement.

### 5. Oscillatory Decay Integral
```
f(x) = exp(-x)·sin(10x)
∫f(x)dx = exp(-x)·(-sin(10x) - 10·cos(10x))/101
```
**Properties:** Oscillatory integrand with decay.

### 6. Near-Singular Power Integral
```
f(x) = x^0.1
∫f(x)dx = x^1.1/1.1
```
**Properties:** Tests handling of near-singular integrands.

## API Reference

### RealFunctionsTestBed Class

```cpp
// Counts
static int getNumFunc();                       // 19 standard functions
static int getNumFuncExtended();               // 8 numerically challenging
static int getNumFuncWithDerivation();         // 11 derivative test functions
static int getNumFuncWithIntegral();           // 16 integration test functions
static int getNumFuncExtendedWithIntegral();   // 6 challenging integrals

// Access by index
const static TestFunctionReal& getFunc(int i);
const static TestFunctionRealWithDerivation& getFuncWithDerivation(int i);
const static TestFunctionRealWithIntegral& getFuncWithIntegral(int i);
const static TestFunctionRealWithDerivation& getFuncExtended(int i);
const static TestFunctionRealWithIntegral& getFuncExtendedWithIntegral(int i);

// Access by name
const static TestFunctionReal& getFunc(const std::string& funcName);
const static TestFunctionRealWithDerivation& getFuncWithDerivation(const std::string& funcName);
const static TestFunctionRealWithIntegral& getFuncWithIntegral(const std::string& funcName);
const static TestFunctionRealWithDerivation& getFuncExtended(const std::string& funcName);
const static TestFunctionRealWithIntegral& getFuncExtendedWithIntegral(const std::string& funcName);

// Category-based access
static std::vector<const TestFunctionReal*> getStandardFunctions();
static std::vector<const TestFunctionRealWithDerivation*> getNumericalAnalysisFunctions();
static std::vector<const TestFunctionRealWithIntegral*> getChallengingIntegrationFunctions();
```

---

# Scalar Functions Test Bed (ℝ³ → ℝ)

## Data Structure

```cpp
template<int N>
struct TestFunctionScalar {
    std::string _funcName;
    ScalarFunction<N> _func;                                    // f(x)
    Real (*_funcDerived)(const VectorN<Real, N>&, int ind);     // ∂f/∂xᵢ
    std::string _funcExpr;
    std::string _funcDerivedExpr;
};
```

## Function Categories

### Quadratic Forms (3 functions)

| Name | f(x,y,z) | ∇f | Properties |
|------|----------|-----|-----------|
| Sphere | x²+y²+z² | (2x, 2y, 2z) | Convex, single minimum at origin |
| Saddle | x²-y²+z² | (2x, -2y, 2z) | Saddle point at origin |
| Weighted Quadratic | x²+4y²+9z² | (2x, 8y, 18z) | Ellipsoidal bowl, tests anisotropic scaling |

### Classic Optimization Functions (3 functions)

| Name | Description | Properties |
|------|-------------|------------|
| **Rosenbrock 3D** | (1-x)²+100(y-x²)²+(1-y)²+100(z-y²)² | Famous "banana" function, narrow curved valley |
| **Rastrigin 3D** | 30 + Σ(xᵢ²-10cos(2πxᵢ)) | Highly multimodal (many local minima) |
| **Ackley 3D** | Complex multimodal | Large nearly flat outer region |

**Purpose:** Standard benchmarks for optimization algorithm testing.

### Physical Potentials (4 functions)

| Name | f(r) | Physical Meaning |
|------|------|------------------|
| **Coulomb Potential** | -1/|r| | Electrostatic/gravitational potential |
| **Harmonic Potential** | ½|r|² | Quantum harmonic oscillator |
| **Lennard-Jones** | 4[(1/r)¹²-(1/r)⁶] | Van der Waals interaction |
| **Morse Potential** | D(1-e^(-a(r-r₀)))² | Molecular bond stretching |

**Purpose:** Tests handling of singularities and physically-motivated functions.

### Additional Mathematical Functions (3 functions)

| Name | f(x,y,z) | ∇f |
|------|----------|-----|
| Gaussian 3D | e^(-(x²+y²+z²)) | -2·(x,y,z)·f |
| Sin Product | sin(x)sin(y)sin(z) | (cos terms with products) |
| Polynomial Cross Terms | x³+y³+z³-3xyz | (3x²-3yz, 3y²-3xz, 3z²-3xy) |

## API Reference

### ScalarFunctionsTestBed Class

```cpp
static int getNumTestFunctionScalar3();  // Returns 15

const static TestFunctionScalar<3>& getTestFunctionScalar3(int i);
const static TestFunctionScalar<3>& getTestFunctionScalar3(const std::string& funcName);
```

---

# Vector Functions Test Bed (ℝ³ → ℝ³)

## Data Structure

```cpp
template<int N>
struct TestFunctionVector {
    std::string _funcName;
    VectorFunction<N> _func;                                         // F(x)
    VectorN<Real, N> (*_funcDerived)(const VectorN<Real, N>&, int);  // Jacobian rows
    std::string _funcExpr;
    std::string _funcDerivedExpr;
    
    // Vector calculus properties
    Real _expectedDiv;       // Expected divergence (NaN if position-dependent)
    Real _expectedCurlMag;   // Expected |curl| (NaN if position-dependent)
    bool _isSolenoidal;      // div F = 0 (incompressible)
    bool _isIrrotational;    // curl F = 0 (conservative)
};
```

## Vector Fields (5 functions)

### 1. Identity Field (Baseline)
```
F(x,y,z) = (x, y, z)
J = I (identity matrix)
div F = 3, curl F = 0
```
**Properties:** Simplest test case, irrotational but not solenoidal.

### 2. Vortex/Rotation Field
```
F(x,y,z) = (-y, x, 0)
J = [[0,-1,0], [1,0,0], [0,0,0]]
div F = 0, curl F = (0, 0, 2)
```
**Properties:** Rigid body rotation around z-axis. Solenoidal (incompressible) with constant vorticity.

### 3. Radial Inverse-Square Field
```
F(x,y,z) = r/|r|³
J_ij = δᵢⱼ/r³ - 3xᵢxⱼ/r⁵
div F = 0 (away from origin), curl F = 0
```
**Properties:** Gravitational/Coulomb field direction. Harmonic field (both solenoidal and irrotational away from singularity).

### 4. Gradient of Product (Conservative)
```
F(x,y,z) = ∇(xyz) = (yz, xz, xy)
J = [[0,z,y], [z,0,x], [y,x,0]] (symmetric)
div F = 0, curl F = 0
```
**Properties:** Conservative field (gradient of scalar). Symmetric Jacobian.

### 5. Complex Mixed Field
```
F(x,y,z) = (x·cos(y)·z², sin(x)·(y²+z²), exp(xy/(z²+1)))
```
**Properties:** Complex function with trigonometric, polynomial, and exponential terms. Position-dependent divergence and curl. Thorough test of numerical differentiation.

## Vector Calculus Properties Summary

| Function | div F | |curl F| | Solenoidal | Irrotational |
|----------|-------|---------|------------|--------------|
| Identity | 3 | 0 | ✗ | ✓ |
| Vortex | 0 | 2 | ✓ | ✗ |
| Radial Inv-Sq | 0 | 0 | ✓ | ✓ |
| Grad(xyz) | 0 | 0 | ✓ | ✓ |
| Complex Mixed | varies | varies | ✗ | ✗ |

## API Reference

### VectorFunctionsTestBed Class

```cpp
static int getNumTestFunctionVector();  // Returns 5

const static TestFunctionVector<3>& getTestFunctionVector(int i);
const static TestFunctionVector<3>& getTestFunctionVector(const std::string& funcName);
```

---

# Usage Examples

## Testing Numerical Differentiation

```cpp
#include "test_data/real_functions_test_bed.h"
#include "MML.h"

using namespace MML;

TEST_CASE("Numerical Derivative Accuracy", "[differentiation]") {
    for (int i = 0; i < TestBeds::RealFunctionsTestBed::getNumFuncWithDerivation(); i++) {
        const auto& f = TestBeds::RealFunctionsTestBed::getFuncWithDerivation(i);
        INFO("Testing derivative of: " << f._funcName);
        
        // Sample point within test interval
        Real x = f._intervalTest->midpoint();
        
        // Use MML's numerical differentiation (4-point central difference)
        Real numDeriv = Derivation::NDer4(f._func, x);
        
        // Compare to analytical derivative from test bed
        Real exactDeriv = f._funcDerived(x);
        
        REQUIRE(std::abs(numDeriv - exactDeriv) < 1e-8);
    }
}

TEST_CASE("Higher Order Derivatives", "[differentiation]") {
    // Test second and third derivatives using complete test functions
    for (int i = 0; i < TestBeds::RealFunctionsTestBed::getNumFunc(); i++) {
        const auto& f = TestBeds::RealFunctionsTestBed::getFunc(i);
        INFO("Testing: " << f._funcName);
        
        Real x = f._intervalTest->midpoint();
        
        // Second derivative
        Real numSecDer = Derivation::NSecDer4(f._func, x);
        Real exactSecDer = f._funcSecDer(x);
        REQUIRE(std::abs(numSecDer - exactSecDer) < 1e-6);
        
        // Third derivative
        Real numThirdDer = Derivation::NThirdDer2(f._func, x);
        Real exactThirdDer = f._funcThirdDer(x);
        REQUIRE(std::abs(numThirdDer - exactThirdDer) < 1e-4);
    }
}
```

## Testing Numerical Integration

```cpp
#include "test_data/real_functions_test_bed.h"
#include "MML.h"

using namespace MML;

TEST_CASE("Integration Accuracy - Simpson", "[integration]") {
    for (int i = 0; i < TestBeds::RealFunctionsTestBed::getNumFuncWithIntegral(); i++) {
        const auto& f = TestBeds::RealFunctionsTestBed::getFuncWithIntegral(i);
        INFO("Integrating: " << f._funcName);
        
        Real a = f._intervalTest->getLowerBound();
        Real b = f._intervalTest->getUpperBound();
        
        // Use MML's Simpson integration
        auto result = IntegrateSimpson(f._func, a, b);
        
        // Compute exact integral from antiderivative
        Real exactIntegral = f._funcIntegrated(b) - f._funcIntegrated(a);
        
        Real relError = std::abs((result.value - exactIntegral) / exactIntegral);
        REQUIRE(relError < 1e-8);
        REQUIRE(result.converged);
    }
}

TEST_CASE("Integration Accuracy - Romberg", "[integration]") {
    for (int i = 0; i < TestBeds::RealFunctionsTestBed::getNumFuncWithIntegral(); i++) {
        const auto& f = TestBeds::RealFunctionsTestBed::getFuncWithIntegral(i);
        
        Real a = f._intervalTest->getLowerBound();
        Real b = f._intervalTest->getUpperBound();
        
        // Romberg integration (higher accuracy)
        auto result = IntegrateRomberg(f._func, a, b);
        Real exactIntegral = f._funcIntegrated(b) - f._funcIntegrated(a);
        
        Real relError = std::abs((result.value - exactIntegral) / exactIntegral);
        REQUIRE(relError < 1e-10);
    }
}

TEST_CASE("Gauss-Legendre Quadrature", "[integration]") {
    for (int i = 0; i < TestBeds::RealFunctionsTestBed::getNumFuncWithIntegral(); i++) {
        const auto& f = TestBeds::RealFunctionsTestBed::getFuncWithIntegral(i);
        
        Real a = f._intervalTest->getLowerBound();
        Real b = f._intervalTest->getUpperBound();
        
        // 10-point Gauss-Legendre quadrature
        auto result = IntegrateGauss10(f._func, a, b);
        Real exactIntegral = f._funcIntegrated(b) - f._funcIntegrated(a);
        
        Real relError = std::abs((result.value - exactIntegral) / exactIntegral);
        REQUIRE(relError < 1e-10);
    }
}
```

## Testing Gradient Computation

```cpp
#include "test_data/scalar_functions_test_bed.h"
#include "MML.h"

using namespace MML;

TEST_CASE("Numerical Gradient Accuracy", "[gradient]") {
    for (int i = 0; i < TestBeds::ScalarFunctionsTestBed::getNumTestFunctionScalar3(); i++) {
        const auto& f = TestBeds::ScalarFunctionsTestBed::getTestFunctionScalar3(i);
        INFO("Testing gradient of: " << f._funcName);
        
        VectorN<Real, 3> point{1.0, 2.0, 3.0};
        
        // Use MML's gradient computation in Cartesian coordinates
        VectorN<Real, 3> numGrad = ScalarFieldOperations::GradientCart<3>(f._func, point);
        
        // Compare to analytical gradient from test bed
        for (int j = 0; j < 3; j++) {
            Real exactPartial = f._funcDerived(point, j);
            REQUIRE(std::abs(numGrad[j] - exactPartial) < 1e-8);
        }
    }
}

TEST_CASE("Gradient with Different Accuracy Orders", "[gradient]") {
    const auto& f = TestBeds::ScalarFunctionsTestBed::getTestFunctionScalar3("Sphere");
    VectorN<Real, 3> point{1.0, 2.0, 3.0};
    
    // Test different derivative orders (2, 4, 6, 8 point formulas)
    auto grad2 = ScalarFieldOperations::GradientCart<3>(f._func, point, 2);
    auto grad4 = ScalarFieldOperations::GradientCart<3>(f._func, point, 4);
    auto grad8 = ScalarFieldOperations::GradientCart<3>(f._func, point, 8);
    
    // Higher order should be more accurate
    VectorN<Real, 3> exact{2.0, 4.0, 6.0};  // Gradient of x²+y²+z² at (1,2,3)
    
    REQUIRE(grad4.IsEqualTo(exact, 1e-10));
    REQUIRE(grad8.IsEqualTo(exact, 1e-12));
}
```

## Testing Vector Field Properties

```cpp
#include "test_data/vector_functions_test_bed.h"
#include "MML.h"

using namespace MML;

TEST_CASE("Divergence Verification", "[vector_calculus]") {
    for (int i = 0; i < TestBeds::VectorFunctionsTestBed::getNumTestFunctionVector(); i++) {
        const auto& F = TestBeds::VectorFunctionsTestBed::getTestFunctionVector(i);
        INFO("Testing divergence of: " << F._funcName);
        
        VectorN<Real, 3> point{0.5, 0.7, 0.9};
        
        // Use MML's divergence computation
        Real numDiv = VectorFieldOperations::DivCart<3>(F._func, point);
        
        // Compare to expected value from test bed (if constant)
        if (!std::isnan(F._expectedDiv)) {
            REQUIRE(std::abs(numDiv - F._expectedDiv) < 1e-8);
        }
    }
}

TEST_CASE("Curl Verification", "[vector_calculus]") {
    for (int i = 0; i < TestBeds::VectorFunctionsTestBed::getNumTestFunctionVector(); i++) {
        const auto& F = TestBeds::VectorFunctionsTestBed::getTestFunctionVector(i);
        INFO("Testing curl of: " << F._funcName);
        
        VectorN<Real, 3> point{0.5, 0.7, 0.9};
        
        // Use MML's curl computation
        VectorN<Real, 3> curl = VectorFieldOperations::CurlCart(F._func, point);
        
        // Check curl magnitude against expected (if constant)
        if (!std::isnan(F._expectedCurlMag)) {
            REQUIRE(std::abs(curl.NormL2() - F._expectedCurlMag) < 1e-8);
        }
        
        // Verify solenoidal/irrotational properties
        if (F._isSolenoidal) {
            Real div = VectorFieldOperations::DivCart<3>(F._func, point);
            REQUIRE(std::abs(div) < 1e-8);
        }
        if (F._isIrrotational) {
            REQUIRE(curl.NormL2() < 1e-8);
        }
    }
}

TEST_CASE("Laplacian of Scalar Fields", "[vector_calculus]") {
    const auto& f = TestBeds::ScalarFunctionsTestBed::getTestFunctionScalar3("Sphere");
    VectorN<Real, 3> point{1.0, 2.0, 3.0};
    
    // Laplacian of x²+y²+z² should be 6 (constant)
    Real lapl = ScalarFieldOperations::LaplacianCart<3>(f._func, point);
    REQUIRE(std::abs(lapl - 6.0) < 1e-8);
}
```

---

# Summary Statistics

| Category | Count | Purpose |
|----------|-------|---------|
| Standard Real Functions | 19 | Complete derivative/integral data |
| Derivative Test Functions | 11 | Focused derivative testing |
| Integration Test Functions | 16 | Focused integration testing |
| Numerically Challenging (Derivatives) | 8 | Stress test differentiation |
| Challenging Integrals | 6 | Stress test integration |
| Scalar Functions (ℝ³→ℝ) | 15 | Gradient and optimization |
| Vector Fields (ℝ³→ℝ³) | 5 | Jacobian, div, curl |
| **Total** | **~80** | |

## References

- Numerical Recipes in C, 3rd Edition - William H. Press et al.
- Test Functions for Optimization - Momin Jamil & Xin-She Yang
- Classical Mathematical Functions - Abramowitz & Stegun
