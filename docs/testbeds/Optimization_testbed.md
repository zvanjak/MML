# Optimization Test Bed Documentation

## Overview

The Optimization Test Bed provides comprehensive test infrastructure for validating MML's optimization algorithms. It includes classic benchmark functions for both 1D and multi-dimensional optimization, covering unimodal, multimodal, and ill-conditioned problems.

**Key Features:**
- **8 1D test functions** (unimodal and multimodal)
- **12+ 2D test functions** (classic benchmarks)
- **N-dimensional test generators** (scalable to arbitrary dimensions)
- Functions with **known global minima** and search domains
- Metadata: category, difficulty, gradient availability, number of minima

---

## Core Files

| File | Purpose |
|------|---------|
| [optimization_test_bed.h](../../test_data/optimization_test_bed.h) | `TestOptimization1D`, `TestOptimizationND<N>` structs, retrieval functions |
| [optimization_defs.h](../../test_data/optimization_defs.h) | All optimization test function implementations |

---

## Part 1: 1D Optimization Test Functions

### Data Structure

```cpp
namespace MML::TestBeds {

/// Test case for 1D optimization algorithms
struct TestOptimization1D {
    std::string name;                        ///< Descriptive name
    std::function<Real(Real)> func;          ///< Objective function f(x)
    std::function<Real(Real)> derivative;    ///< f'(x) for gradient-based methods
    Real trueMinimumX;                       ///< Known x* location of minimum
    Real trueMinimumF;                       ///< Known f(x*) value at minimum
    Real searchLow;                          ///< Lower bound of search range
    Real searchHigh;                         ///< Upper bound of search range
    
    // Metadata
    std::string category;                    ///< unimodal, multimodal, challenging
    std::string description;                 ///< Detailed description
    bool hasDerivative = true;               ///< Whether derivative is available
    bool isMultimodal = false;               ///< Has multiple local minima
    int difficulty = 1;                      ///< 1=easy, 2=medium, 3=hard
};

}
```

### 1D Test Functions Catalog (8 total)

#### Unimodal Functions (5)

| # | Name | Function f(x) | Domain | x* | f(x*) | Difficulty |
|---|------|---------------|--------|----|----|------------|
| 1 | **Quadratic** | (x-2)² + 1 | [-5, 10] | 2.0 | 1.0 | Easy |
| 2 | **Quartic** | (x-1)⁴ | [-3, 5] | 1.0 | 0.0 | Medium (flat) |
| 3 | **Cosh** | exp(x-3) + exp(-(x-3)) | [0, 6] | 3.0 | 2.0 | Easy |
| 4 | **AbsValue** | |x-1.5| + 2 | [-2, 5] | 1.5 | 2.0 | Medium (non-smooth) |
| 5 | **LogBarrier** | x² - ln(x) | [0.1, 5] | 1/√2 | 0.847 | Medium |

#### Multimodal Functions (3)

| # | Name | Function f(x) | Domain | x* | f(x*) | Local Minima |
|---|------|---------------|--------|----|----|--------------|
| 6 | **Rastrigin 1D** | 10 + x² - 10cos(2πx) | [-5.12, 5.12] | 0.0 | 0.0 | Many |
| 7 | **Ackley 1D** | -20exp(-0.2|x|) - exp(cos(2πx)) + 20 + e | [-32.768, 32.768] | 0.0 | 0.0 | Many |
| 8 | **Sinusoidal** | sin(x) + sin(3x)/3 + 0.5 | [-π, 2π] | ≈-1.9 | ≈-0.37 | Several |

---

## Part 2: N-Dimensional Optimization Test Functions

### Data Structure

```cpp
namespace MML::TestBeds {

/// Test case for N-dimensional optimization algorithms
template<int N>
struct TestOptimizationND {
    std::string name;                                      ///< Descriptive name
    std::function<Real(const VectorN<Real, N>&)> func;     ///< Objective function
    std::function<VectorN<Real, N>(const VectorN<Real, N>&)> gradient;  ///< ∇f
    std::vector<VectorN<Real, N>> trueMinima;              ///< Known minimum locations
    Real trueMinimumF;                                     ///< Known minimum value
    VectorN<Real, N> domainLow;                            ///< Lower bounds
    VectorN<Real, N> domainHigh;                           ///< Upper bounds
    VectorN<Real, N> suggestedStart;                       ///< Suggested starting point
    
    // Metadata
    std::string category;                                  ///< convex, valley, multimodal, etc.
    std::string description;                               ///< Detailed description
    bool hasGradient = false;                              ///< Whether gradient is available
    bool isMultimodal = false;                             ///< Has multiple local minima
    int numMinima = 1;                                     ///< Number of global minima
    int difficulty = 1;                                    ///< 1=easy, 2=medium, 3=hard, 4=very hard
};

}
```

### 2D Test Functions Catalog (12 total)

#### Convex/Unimodal Functions

| # | Name | Category | x* | f(x*) | Difficulty | Description |
|---|------|----------|----|----|------------|-------------|
| 1 | **Sphere 2D** | Convex | (0, 0) | 0 | Easy | Simplest benchmark |
| 2 | **Booth** | Convex | (1, 3) | 0 | Easy | Quadratic bowl |
| 3 | **Matyas** | Convex | (0, 0) | 0 | Easy | Plate-shaped |
| 4 | **Beale** | Unimodal | (3, 0.5) | 0 | Medium | Flat regions |
| 5 | **Easom** | Unimodal | (π, π) | -1 | Hard | Nearly flat with sharp minimum |

#### Valley Functions

| # | Name | Category | x* | f(x*) | Difficulty | Description |
|---|------|----------|----|----|------------|-------------|
| 6 | **Rosenbrock** | Valley | (1, 1) | 0 | Medium | Classic banana-shaped valley |

#### Multimodal Functions

| # | Name | # Minima | Global Minimum | f(x*) | Difficulty |
|---|------|----------|---------------|-------|------------|
| 7 | **Himmelblau** | 4 | (3,2), (-2.8,3.1), (-3.8,-3.3), (3.6,-1.8) | 0 | Medium |
| 8 | **Goldstein-Price** | 1 (several local) | (0, -1) | 3 | Hard |
| 9 | **Three-Hump Camel** | 1 | (0, 0) | 0 | Medium |
| 10 | **Six-Hump Camel** | 2 | (±0.09, ∓0.71) | -1.032 | Medium |
| 11 | **Rastrigin 2D** | 1 (~100 local) | (0, 0) | 0 | Hard |
| 12 | **Ackley 2D** | 1 (many local) | (0, 0) | 0 | Hard |

### Higher-Dimensional Test Functions

| Function | Dim | x* | f(x*) | Has Gradient | Category |
|----------|-----|----|----|--------------|----------|
| **Sphere ND** | N | (0,...,0) | 0 | ✅ | Convex |
| **Rosenbrock ND** | N | (1,...,1) | 0 | ✅ | Valley |
| **Rastrigin ND** | N | (0,...,0) | 0 | ❌ | Multimodal |
| **Ackley ND** | N | (0,...,0) | 0 | ❌ | Multimodal |
| **Powell 4D** | 4 | (0,0,0,0) | 0 | ❌ | Ill-conditioned |

### Additional Functions in Definitions

| # | Function | Formula | x* | Category |
|---|----------|---------|-----|----------|
| 13 | Sum of Squares | Σ i·xᵢ² | 0 | Convex |
| 14 | Rotated Hyper-Ellipsoid | Σᵢ(Σⱼ≤ᵢ xⱼ)² | 0 | Convex (non-separable) |
| 15 | Trid | Σ(xᵢ-1)² - Σxᵢxᵢ₋₁ | varies | Convex |
| 16 | Griewank | 1 + Σxᵢ²/4000 - Πcos(xᵢ/√i) | 0 | Multimodal |
| 17 | Schwefel | 418.98n - Σxᵢsin(√|xᵢ|) | 420.97 | Deceptive |
| 18 | Levy | Complex | 1 | Multimodal |
| 19 | Styblinski-Tang | 0.5Σ(xᵢ⁴-16xᵢ²+5xᵢ) | -2.9 | Multimodal |
| 20 | Dixon-Price | (x₁-1)² + Σi(2xᵢ²-xᵢ₋₁)² | varies | Ill-conditioned |
| 21 | Zakharov | Σxᵢ² + (Σ0.5i·xᵢ)² + (...)⁴ | 0 | Ill-conditioned |

---

## Search Domain Bounds

Standard search domains for benchmark functions:

| Function | Low | High |
|----------|-----|------|
| Sphere | -5.12 | 5.12 |
| Rosenbrock | -5.0 | 10.0 |
| Rastrigin | -5.12 | 5.12 |
| Ackley | -32.768 | 32.768 |
| Griewank | -600.0 | 600.0 |
| Schwefel | -500.0 | 500.0 |
| Beale | -4.5 | 4.5 |
| Goldstein-Price | -2.0 | 2.0 |
| Himmelblau | -5.0 | 5.0 |

---

## API Reference

### Retrieval Functions - 1D

```cpp
namespace MML::TestBeds {

// Individual test getters
TestOptimization1D getQuadratic1DTest();
TestOptimization1D getQuartic1DTest();
TestOptimization1D getCosh1DTest();
TestOptimization1D getAbsValue1DTest();
TestOptimization1D getLogBarrier1DTest();
TestOptimization1D getRastrigin1DTest();
TestOptimization1D getAckley1DTest();
TestOptimization1D getSinusoidal1DTest();

// Collection getters
std::vector<TestOptimization1D> getAll1DOptimizationTests();     // All 8 tests
std::vector<TestOptimization1D> getUnimodal1DTests();            // 5 unimodal
std::vector<TestOptimization1D> getMultimodal1DTests();          // 3 multimodal

}
```

### Retrieval Functions - 2D

```cpp
namespace MML::TestBeds {

// Individual test getters
TestOptimizationND<2> getSphere2DTest();
TestOptimizationND<2> getRosenbrock2DTest();
TestOptimizationND<2> getBeale2DTest();
TestOptimizationND<2> getBooth2DTest();
TestOptimizationND<2> getMatyas2DTest();
TestOptimizationND<2> getHimmelblau2DTest();
TestOptimizationND<2> getGoldsteinPrice2DTest();
TestOptimizationND<2> getThreeHumpCamel2DTest();
TestOptimizationND<2> getSixHumpCamel2DTest();
TestOptimizationND<2> getEasom2DTest();
TestOptimizationND<2> getRastrigin2DTest();
TestOptimizationND<2> getAckley2DTest();

// Collection getters
std::vector<TestOptimizationND<2>> getAll2DOptimizationTests();       // All 12
std::vector<TestOptimizationND<2>> getConvex2DTests();                // Sphere, Booth, Matyas
std::vector<TestOptimizationND<2>> getMultimodal2DTests();            // Himmelblau, etc.
std::vector<TestOptimizationND<2>> getClassicBenchmark2DTests();      // Famous ones

}
```

### Retrieval Functions - N-Dimensional

```cpp
namespace MML::TestBeds {

// Template generators (any dimension)
template<int N> TestOptimizationND<N> getSphereNDTest();
template<int N> TestOptimizationND<N> getRosenbrockNDTest();
template<int N> TestOptimizationND<N> getRastriginNDTest();
template<int N> TestOptimizationND<N> getAckleyNDTest();

// Specific higher-dimensional
TestOptimizationND<4> getPowell4DTest();

}
```

### Verification Utilities

```cpp
namespace MML::TestBeds {

/// Verify a 1D result is close to true minimum
bool verify1DMinimum(const TestOptimization1D& test, Real foundX, Real tolerance = 1e-6);

/// Verify an ND result (checks if close to ANY known minimum)
template<int N>
bool verifyNDMinimum(const TestOptimizationND<N>& test, const VectorN<Real, N>& found, 
                     Real tolerance = 1e-6);

/// Compute distance to nearest known minimum
template<int N>
Real computeMinimumError(const TestOptimizationND<N>& test, const VectorN<Real, N>& found);

/// Create custom test case
TestOptimization1D createCustom1DTest(const std::string& name, std::function<Real(Real)> func,
                                       Real trueMinX, Real trueMinF, Real searchLow, Real searchHigh,
                                       const std::string& description = "Custom test case");

/// Generate random starting point within domain
template<int N>
VectorN<Real, N> generateRandomStart(const TestOptimizationND<N>& test, unsigned int seed = 42);

}
```

### Function Wrappers

```cpp
namespace MML::TestBeds {

/// Adapt std::function to IRealFunction (1D)
class OptRealFunctionWrapper : public IRealFunction {
public:
    OptRealFunctionWrapper(std::function<Real(Real)> f);
    Real operator()(Real x) const override;
};

/// Adapt std::function to IScalarFunction<N> (ND)
template<int N>
class OptScalarFunctionWrapper : public IScalarFunction<N> {
public:
    OptScalarFunctionWrapper(std::function<Real(const VectorN<Real, N>&)> f);
    Real operator()(const VectorN<Real, N>& x) const override;
};

}
```

---

## Usage Examples

### Example 1: Test Golden Section Search (1D)

```cpp
#include "test_data/optimization_test_bed.h"
#include "algorithms/Optimization.h"

using namespace MML;
using namespace MML::TestBeds;

void testGoldenSection() {
    // Get all unimodal 1D tests
    auto tests = getUnimodal1DTests();
    
    for (const auto& test : tests) {
        // Wrap function for MML interface
        OptRealFunctionWrapper func(test.func);
        
        // Bracket the minimum
        auto bracket = Minimization::BracketMinimum(func, test.searchLow, test.searchHigh);
        
        // Find minimum using golden section
        auto result = Minimization::GoldenSectionSearch(func, bracket, 1e-8);
        
        // Verify result
        bool success = verify1DMinimum(test, result.xmin, 1e-6);
        std::cout << test.name << ": " << (success ? "PASS" : "FAIL")
                  << " (error = " << std::abs(result.xmin - test.trueMinimumX) << ")" << std::endl;
    }
}
```

### Example 2: Test Brent's Method with Derivatives (1D)

```cpp
void testBrentWithDerivative() {
    auto tests = getAll1DOptimizationTests();
    
    for (const auto& test : tests) {
        if (!test.hasDerivative) continue;  // Skip non-differentiable
        
        OptRealFunctionWrapper func(test.func);
        RealFunctionFromStdFunc deriv(test.derivative);
        
        auto bracket = Minimization::BracketMinimum(func, test.searchLow, test.searchHigh);
        auto result = Minimization::BrentMinimizeWithDeriv(func, deriv, bracket, 1e-10);
        
        Real error = std::abs(result.xmin - test.trueMinimumX);
        std::cout << test.name << ": x* = " << result.xmin 
                  << ", error = " << error 
                  << ", iters = " << result.iterations << std::endl;
    }
}
```

### Example 3: Test Nelder-Mead on 2D Functions

```cpp
#include "test_data/optimization_test_bed.h"
#include "algorithms/Optimization/OptimizationMultidim.h"

void testNelderMead2D() {
    auto tests = getAll2DOptimizationTests();
    
    for (const auto& test : tests) {
        // Wrap function
        OptScalarFunctionWrapper<2> func(test.func);
        
        // Run Nelder-Mead from suggested starting point
        auto result = NelderMeadMinimize(func, test.suggestedStart, 1.0, 1e-8);
        
        // Check if found any global minimum
        bool success = verifyNDMinimum(test, result.minLoc, 1e-4);
        Real error = computeMinimumError(test, result.minLoc);
        
        std::cout << test.name 
                  << ": " << (success ? "PASS" : "FAIL")
                  << ", f(x*) = " << result.minVal
                  << ", error = " << error
                  << ", iters = " << result.iterations << std::endl;
    }
}
```

### Example 4: Test Powell's Method on Rosenbrock

```cpp
void testPowellRosenbrock() {
    auto test = getRosenbrock2DTest();
    OptScalarFunctionWrapper<2> func(test.func);
    
    // Powell's method (derivative-free)
    Powell optimizer(1e-10, 500);
    auto result = optimizer.minimize(func, test.suggestedStart);
    
    std::cout << "Rosenbrock via Powell:" << std::endl;
    std::cout << "  x* = (" << result.minLoc[0] << ", " << result.minLoc[1] << ")" << std::endl;
    std::cout << "  f(x*) = " << result.minVal << std::endl;
    std::cout << "  Expected: (1, 1), f = 0" << std::endl;
    std::cout << "  Error: " << computeMinimumError(test, result.minLoc) << std::endl;
}
```

### Example 5: Test Conjugate Gradient with Gradient

```cpp
void testConjugateGradient() {
    // Get tests with gradients
    auto tests = getAll2DOptimizationTests();
    
    for (const auto& test : tests) {
        if (!test.hasGradient) continue;
        
        OptScalarFunctionWrapper<2> func(test.func);
        
        // Create gradient wrapper
        auto gradFunc = [&test](const VectorN<Real, 2>& x) -> VectorN<Real, 2> {
            return test.gradient(x);
        };
        
        ConjugateGradient cg(1e-10, 1e-10, 500);
        auto result = cg.minimize(func, gradFunc, test.suggestedStart);
        
        bool success = verifyNDMinimum(test, result.minLoc, 1e-5);
        std::cout << test.name << " (CG): " << (success ? "PASS" : "FAIL") << std::endl;
    }
}
```

### Example 6: Algorithm Comparison on Classic Benchmarks

```cpp
void compareAlgorithms() {
    auto benchmarks = getClassicBenchmark2DTests();
    
    std::cout << "Algorithm Comparison:\n";
    std::cout << std::setw(20) << "Function" 
              << std::setw(15) << "Nelder-Mead"
              << std::setw(15) << "Powell" << std::endl;
    
    for (const auto& test : benchmarks) {
        OptScalarFunctionWrapper<2> func(test.func);
        
        // Nelder-Mead
        NelderMead nm(1e-8, 5000);
        auto nmResult = nm.minimize(func, test.suggestedStart);
        Real nmError = computeMinimumError(test, nmResult.minLoc);
        
        // Powell
        Powell pw(1e-8, 500);
        auto pwResult = pw.minimize(func, test.suggestedStart);
        Real pwError = computeMinimumError(test, pwResult.minLoc);
        
        std::cout << std::setw(20) << test.name
                  << std::setw(15) << nmError
                  << std::setw(15) << pwError << std::endl;
    }
}
```

### Example 7: Scale to Higher Dimensions

```cpp
void testHighDimensional() {
    // Test in 5D, 10D, 20D
    for (int N : {5, 10, 20}) {
        if (N == 5) {
            auto test = getSphereNDTest<5>();
            std::cout << "Testing " << test.name << std::endl;
        } else if (N == 10) {
            auto test = getRosenbrockNDTest<10>();
            std::cout << "Testing " << test.name << " (difficulty: " << test.difficulty << ")" << std::endl;
        } else {
            auto test = getRastriginNDTest<20>();
            std::cout << "Testing " << test.name << " (multimodal, ~10^20 local minima)" << std::endl;
        }
    }
}
```

---

## MML Optimization Solver Reference

### 1D Optimization Methods

| Method | Derivative? | Order | Use Case |
|--------|------------|-------|----------|
| `Minimization::BracketMinimum()` | No | - | Find initial bracket |
| `Minimization::GoldenSectionSearch()` | No | Linear | Robust, non-smooth |
| `Minimization::BrentMinimize()` | No | Superlinear | **Recommended default** |
| `Minimization::BrentMinimizeWithDeriv()` | Yes | Superlinear | Faster if gradient available |

### N-Dimensional Optimization Methods

| Method | Gradient? | Order | Use Case |
|--------|----------|-------|----------|
| `NelderMead` | No | Sublinear | Robust, general purpose |
| `Powell` | No | Superlinear | Derivative-free, valleys |
| `ConjugateGradient` | Yes | Superlinear | Large-scale, smooth |
| `BFGS` | Yes | Superlinear | **Recommended if gradient available** |

### Algorithm Selection Guide

```
Is gradient available?
├── No → Use derivative-free method
│   ├── 1D problem? → BrentMinimize
│   ├── 2D-10D, simple? → NelderMead
│   └── Valley structure? → Powell
│
└── Yes → Use gradient-based method
    ├── 1D? → BrentMinimizeWithDeriv
    ├── Small N? → BFGS
    └── Large N (>100)? → ConjugateGradient
```

---

## Test Function Categories

### By Difficulty

| Difficulty | 1D Functions | 2D Functions |
|------------|--------------|--------------|
| **Easy (1)** | Quadratic, Cosh | Sphere, Booth, Matyas |
| **Medium (2)** | Quartic, AbsValue, LogBarrier, Sinusoidal | Rosenbrock, Beale, Himmelblau, Camel |
| **Hard (3)** | Rastrigin, Ackley | Goldstein-Price, Easom, Rastrigin, Ackley |
| **Very Hard (4)** | - | High-dim Rastrigin |

### By Category

| Category | Description | Examples |
|----------|-------------|----------|
| **Convex** | Single minimum, no local minima | Sphere, Booth, Matyas |
| **Valley** | Long curved valley | Rosenbrock |
| **Unimodal** | Single minimum but complex shape | Beale, Easom |
| **Multimodal** | Multiple local minima | Rastrigin, Ackley, Himmelblau |
| **Ill-conditioned** | Narrow valleys, poor scaling | Powell, Dixon-Price |
| **Deceptive** | Global minimum far from local | Schwefel |

---

## See Also

- [RootFinding_testbed.md](RootFinding_testbed.md) - Root finding test cases
- [Functions_testbed.md](Functions_testbed.md) - General test functions
- [Optimization.h](../../mml/algorithms/Optimization.h) - 1D optimization algorithms
- [OptimizationMultidim.h](../../mml/algorithms/Optimization/OptimizationMultidim.h) - N-D algorithms