# Root Finding Test Bed Documentation

## Overview

The Root Finding Test Bed provides comprehensive test infrastructure for validating MML's root finding algorithms. It includes polynomial equations, transcendental equations, multiple root cases, and challenging pathological functions with known exact roots.

**Key Features:**
- **25 test functions** across 7 categories
- Functions with **known analytical roots** and bracketing intervals
- Metadata: category, difficulty, multiplicity, singularity information
- **High-precision root constants** for verification

---

## Core Files

| File | Purpose |
|------|---------|
| [root_finding_test_bed.h](../../test_data/root_finding_test_bed.h) | `TestRootFunction` struct, retrieval functions, verification utilities |
| [root_finding_defs.h](../../test_data/root_finding_defs.h) | All test function implementations and root constants |

---

## Data Structure

```cpp
namespace MML::TestBeds {

/// Test case for root finding algorithms
struct TestRootFunction {
    std::string name;                              ///< Descriptive name
    std::function<Real(Real)> func;                ///< The function f(x)
    std::function<Real(Real)> derivative;          ///< f'(x) for Newton's method
    std::vector<Real> knownRoots;                  ///< Analytically known root values
    std::vector<int> multiplicities;               ///< Multiplicity of each root (1=simple, 2=double, etc.)
    std::vector<std::pair<Real,Real>> brackets;    ///< [low, high] bracketing intervals
    
    // Metadata
    std::string category;                          ///< polynomial, transcendental, multiple_roots, etc.
    std::string description;                       ///< Detailed description
    bool hasDerivative = true;                     ///< Whether derivative is available
    bool hasSingularity = false;                   ///< Whether function has singularities
    Real singularityLocation = 0;                  ///< Location of singularity if present
    int difficulty = 1;                            ///< 1=easy, 2=medium, 3=hard, 4=challenging
};

}
```

---

## Test Functions Catalog (25 total)

### Category 1: Polynomial Functions (5)

Well-conditioned polynomials with exact roots.

| # | Name | Function f(x) | Roots | Difficulty |
|---|------|---------------|-------|------------|
| 1 | **Quadratic √2** | x² - 2 | ±√2 ≈ ±1.4142 | Easy |
| 2 | **Cubic 1,2,3** | x³ - 6x² + 11x - 6 | 1, 2, 3 | Easy |
| 3 | **Quadratic 2,3** | x² - 5x + 6 | 2, 3 | Easy |
| 4 | **Quartic ±1,±2** | x⁴ - 5x² + 4 | -2, -1, 1, 2 | Medium |
| 5 | **Wilkinson-5** | (x-1)(x-2)(x-3)(x-4)(x-5) | 1, 2, 3, 4, 5 | Medium |

### Category 2: Transcendental Equations (6)

Classical benchmarks with mathematical constants.

| # | Name | Function f(x) | Root | Significance |
|---|------|---------------|------|--------------|
| 6 | **Dottie Number** | x - cos(x) | 0.7391 | Mathematical constant |
| 7 | **Kepler Equation** | x - 0.5sin(x) - 0.5 | 0.7582 | Orbital mechanics |
| 8 | **Omega Constant** | exp(-x) - x | 0.5671 | Lambert W related |
| 9 | **Euler Number** | ln(x) - 1 | e ≈ 2.7183 | Singularity at x=0 |
| 10 | **Sin Half** | sin(x) - 0.5 | π/6 ≈ 0.5236 | Trigonometric |
| 11 | **Tangent Intersection** | tan(x) - x | 0, 4.4934, ... | Multiple roots, singularities |

### Category 3: Multiple/Repeated Roots (3)

Challenging for most methods due to reduced convergence order.

| # | Name | Function f(x) | Roots | Multiplicity |
|---|------|---------------|-------|--------------|
| 12 | **Double Root** | (x-1)² | 1 | 2 |
| 13 | **Triple Root** | (x-2)³ | 2 | 3 |
| 14 | **Mixed Multiplicity** | (x-1)²(x-3) | 1 (×2), 3 (×1) | Mixed |

### Category 4: Closely Spaced Roots (2)

Ill-conditioned problems testing precision.

| # | Name | Function f(x) | Roots | Separation |
|---|------|---------------|-------|------------|
| 15 | **Close Roots 1** | (x-1)(x-1.001) | 1.0, 1.001 | 0.001 |
| 16 | **Close Roots 2** | (x-2)(x-2.0001) | 2.0, 2.0001 | 0.0001 |

### Category 5: Pathological/Challenging Cases (4)

Stress tests for algorithm robustness.

| # | Name | Function f(x) | Root | Challenge |
|---|------|---------------|------|-----------|
| 17 | **Steep Exponential** | exp(x) - 10000 | ln(10000) ≈ 9.21 | Large derivative |
| 18 | **Flat Cubic** | (x-1)³ - 0.001 | 1 + ∛0.001 ≈ 1.1 | Very flat near root |
| 19 | **Oscillatory Decay** | sin(10x)·exp(-x) | kπ/10 for k=0,1,2,... | Many roots |
| 20 | **Steep Atan** | atan(1000(x-1)) | 1.0 | Derivative ~1000 |

### Category 6: Functions with Singularities (2)

Roots near singularities require careful handling.

| # | Name | Function f(x) | Root | Singularity |
|---|------|---------------|------|-------------|
| 21 | **1/x - 2** | 1/x - 2 | 0.5 | x = 0 |
| 22 | **Log Singularity** | ln(x) + x - 2 | ≈1.557 | x = 0 |

### Category 7: Physics/Engineering Applications (3)

Real-world problems from science and engineering.

| # | Name | Equation | Root | Application |
|---|------|----------|------|-------------|
| 23 | **Van der Waals** | V³ - 2V² + V - 0.1 | ≈0.109 | Real gas thermodynamics |
| 24 | **Planck Radiation** | x - 5(1 - exp(-x)) | ≈4.965 | Wien displacement law |
| 25 | **Lambert W** | x·exp(x) - 1 | W(1) ≈ 0.5671 | Combinatorics, delay equations |

---

## High-Precision Root Constants

```cpp
namespace MML::TestBeds {

// Mathematical constants (50+ digit precision)
const Real ROOT_SQRT2     = 1.41421356237309504880168872420969807856967187537694;
const Real ROOT_DOTTIE    = 0.73908513321516064165531208767387340401341175890076;
const Real ROOT_OMEGA     = 0.56714329040978387299996866221035554975381578718651;
const Real ROOT_E         = 2.71828182845904523536028747135266249775724709369995;
const Real ROOT_PI_6      = 0.52359877559829887307710723054658381403286156656252;
const Real ROOT_TAN_1     = 4.49340945790906417535833613264473925167866328789787;
const Real ROOT_LAMBERT_W1= 0.56714329040978387299996866221035554975381578718651;
const Real ROOT_WIEN      = 4.96511423174427630226686427576612621889222458284869;
const Real ROOT_KEPLER    = 0.75818471929987481918927862261787169523315279313234;
const Real ROOT_LN10000   = 9.21034037197618203544374908830428643115659825936887;

}
```

---

## API Reference

### Retrieval Functions

```cpp
namespace MML::TestBeds {

// Get all tests
std::vector<TestRootFunction> getAllRootFindingTests();         // All 25 tests

// By category
std::vector<TestRootFunction> getPolynomialRootTests();         // 5 polynomial
std::vector<TestRootFunction> getTranscendentalRootTests();     // 6 transcendental
std::vector<TestRootFunction> getMultipleRootTests();           // 3 multiple roots
std::vector<TestRootFunction> getCloseRootTests();              // 2 close roots
std::vector<TestRootFunction> getPathologicalRootTests();       // 4 pathological
std::vector<TestRootFunction> getSingularityRootTests();        // 2 singularity
std::vector<TestRootFunction> getPhysicsRootTests();            // 3 physics/engineering

// By difficulty
std::vector<TestRootFunction> getEasyRootTests();               // Difficulty 1
std::vector<TestRootFunction> getChallengingRootTests();        // Difficulty 3+

}
```

### Individual Test Getters

```cpp
namespace MML::TestBeds {

// Polynomials
TestRootFunction getQuadraticSqrt2Test();
TestRootFunction getCubic123Test();
TestRootFunction getQuadratic23Test();
TestRootFunction getQuartic1122Test();
TestRootFunction getWilkinson5Test();

// Transcendental
TestRootFunction getDottieNumberTest();
TestRootFunction getKeplerEquationTest();
TestRootFunction getOmegaConstantTest();
TestRootFunction getEulerNumberTest();
TestRootFunction getSinHalfTest();
TestRootFunction getTangentIntersectionTest();

// Multiple roots
TestRootFunction getDoubleRootTest();
TestRootFunction getTripleRootTest();
TestRootFunction getMixedMultiplicityTest();

// Close roots
TestRootFunction getCloseRoots1Test();
TestRootFunction getCloseRoots2Test();

// Pathological
TestRootFunction getSteepExponentialTest();
TestRootFunction getFlatCubicTest();
TestRootFunction getOscillatoryDecayTest();
TestRootFunction getSteepAtanTest();

// Singularity
TestRootFunction getSingularityNearRootTest();
TestRootFunction getLogSingularityTest();

// Physics
TestRootFunction getVanDerWaalsTest();
TestRootFunction getPlanckRadiationTest();
TestRootFunction getLambertWTest();

}
```

### Verification Utilities

```cpp
namespace MML::TestBeds {

/// Verify a computed root against test case
bool verifyRoot(const TestRootFunction& test, Real computedRoot, 
                size_t rootIndex = 0, Real tolerance = 1e-10);

/// Compute error of root approximation
Real computeRootError(const TestRootFunction& test, Real computedRoot, 
                      size_t rootIndex = 0);

/// Evaluate |f(x)| at computed root (should be ~0)
Real evaluateFunctionAtRoot(const TestRootFunction& test, Real computedRoot);

/// Create custom test case
TestRootFunction createCustomRootTest(const std::string& name,
    std::function<Real(Real)> func, std::function<Real(Real)> derivative,
    Real knownRoot, Real bracketLow, Real bracketHigh,
    const std::string& description = "Custom test case");

}
```

### Function Wrapper

```cpp
namespace MML::TestBeds {

/// Adapt std::function to IRealFunction interface
class RealFunctionWrapper : public IRealFunction {
public:
    RealFunctionWrapper(std::function<Real(Real)> f);
    Real operator()(Real x) const override;
};

}
```

---

## Usage Examples

### Example 1: Test Bisection Method

```cpp
#include "test_data/root_finding_test_bed.h"
#include "algorithms/RootFinding.h"

using namespace MML;
using namespace MML::TestBeds;

void testBisection() {
    auto tests = getPolynomialRootTests();
    
    for (const auto& test : tests) {
        RealFunctionWrapper func(test.func);
        
        // Test each known root
        for (size_t i = 0; i < test.knownRoots.size(); ++i) {
            Real low = test.brackets[i].first;
            Real high = test.brackets[i].second;
            
            Real root = RootFinding::FindRootBisection(func, low, high, 1e-10);
            
            Real error = computeRootError(test, root, i);
            std::cout << test.name << " root[" << i << "]: " 
                      << root << ", error = " << error << std::endl;
        }
    }
}
```

### Example 2: Test Newton-Raphson

```cpp
void testNewtonRaphson() {
    auto tests = getAllRootFindingTests();
    
    for (const auto& test : tests) {
        if (!test.hasDerivative) continue;  // Skip if no derivative
        if (test.difficulty > 2) continue;  // Skip hard cases initially
        
        RealFunctionWrapper func(test.func);
        
        for (size_t i = 0; i < test.knownRoots.size(); ++i) {
            Real low = test.brackets[i].first;
            Real high = test.brackets[i].second;
            
            try {
                Real root = RootFinding::FindRootNewton(func, low, high, 1e-12);
                bool success = verifyRoot(test, root, i, 1e-10);
                
                std::cout << test.name << ": " << (success ? "PASS" : "FAIL")
                          << ", root = " << root << std::endl;
            } catch (const RootFindingError& e) {
                std::cout << test.name << ": EXCEPTION - " << e.what() << std::endl;
            }
        }
    }
}
```

### Example 3: Test Brent's Method (Gold Standard)

```cpp
void testBrent() {
    // Brent's method works well on challenging cases
    auto tests = getChallengingRootTests();
    
    for (const auto& test : tests) {
        RealFunctionWrapper func(test.func);
        
        for (size_t i = 0; i < test.knownRoots.size(); ++i) {
            Real low = test.brackets[i].first;
            Real high = test.brackets[i].second;
            
            Real root = RootFinding::FindRootBrent(func, low, high, 1e-12);
            Real fval = evaluateFunctionAtRoot(test, root);
            
            std::cout << test.name << ": root = " << root 
                      << ", |f(root)| = " << fval << std::endl;
        }
    }
}
```

### Example 4: Compare All Methods

```cpp
void compareAllMethods() {
    auto test = getDottieNumberTest();  // Classic benchmark
    RealFunctionWrapper func(test.func);
    Real low = test.brackets[0].first;
    Real high = test.brackets[0].second;
    
    std::cout << "Dottie Number Test (true root = " << ROOT_DOTTIE << ")\n";
    std::cout << std::setw(20) << "Method" 
              << std::setw(20) << "Root"
              << std::setw(15) << "Error" << std::endl;
    
    Real root_bisec = RootFinding::FindRootBisection(func, low, high, 1e-12);
    Real root_newton = RootFinding::FindRootNewton(func, low, high, 1e-12);
    Real root_secant = RootFinding::FindRootSecant(func, low, high, 1e-12);
    Real root_ridders = RootFinding::FindRootRidders(func, low, high, 1e-12);
    Real root_brent = RootFinding::FindRootBrent(func, low, high, 1e-12);
    Real root_falsepos = RootFinding::FindRootFalsePosition(func, low, high, 1e-12);
    
    auto printResult = [&](const char* name, Real root) {
        std::cout << std::setw(20) << name 
                  << std::setw(20) << std::setprecision(15) << root
                  << std::setw(15) << std::abs(root - ROOT_DOTTIE) << std::endl;
    };
    
    printResult("Bisection", root_bisec);
    printResult("Newton-Raphson", root_newton);
    printResult("Secant", root_secant);
    printResult("Ridders", root_ridders);
    printResult("Brent", root_brent);
    printResult("False Position", root_falsepos);
}
```

### Example 5: Test Multiple Roots

```cpp
void testMultipleRoots() {
    auto tests = getMultipleRootTests();
    
    std::cout << "Multiple Root Tests (Newton convergence is slower)\n";
    
    for (const auto& test : tests) {
        RealFunctionWrapper func(test.func);
        Real low = test.brackets[0].first;
        Real high = test.brackets[0].second;
        
        std::cout << "\n" << test.name 
                  << " (multiplicity " << test.multiplicities[0] << "):\n";
        
        // Bisection: always works, linear convergence
        Real root_bisec = RootFinding::FindRootBisection(func, low, high, 1e-10);
        
        // Brent: robust, still converges well
        Real root_brent = RootFinding::FindRootBrent(func, low, high, 1e-10);
        
        std::cout << "  Bisection: " << root_bisec 
                  << ", error = " << std::abs(root_bisec - test.knownRoots[0]) << std::endl;
        std::cout << "  Brent:     " << root_brent 
                  << ", error = " << std::abs(root_brent - test.knownRoots[0]) << std::endl;
    }
}
```

### Example 6: Find All Roots of a Polynomial

```cpp
void findAllPolynomialRoots() {
    auto test = getWilkinson5Test();  // 5 roots: 1, 2, 3, 4, 5
    RealFunctionWrapper func(test.func);
    
    // Use bracket finding to locate all roots
    Vector<Real> bracketLow, bracketHigh;
    int numRoots = RootFinding::FindRootBrackets(func, 0.0, 6.0, 100, 
                                                  bracketLow, bracketHigh);
    
    std::cout << "Wilkinson-5 polynomial: found " << numRoots << " brackets\n";
    
    for (int i = 0; i < numRoots; ++i) {
        Real root = RootFinding::FindRootBrent(func, bracketLow[i], bracketHigh[i], 1e-12);
        std::cout << "  Root " << i << ": " << root << std::endl;
    }
}
```

### Example 7: Physics Application - Kepler Equation

```cpp
void solveKeplerEquation() {
    // Kepler's equation: E - e*sin(E) = M
    // Given: eccentricity e = 0.5, mean anomaly M = 0.5
    // Find: eccentric anomaly E
    
    auto test = getKeplerEquationTest();
    RealFunctionWrapper func(test.func);
    
    Real low = test.brackets[0].first;
    Real high = test.brackets[0].second;
    
    // Brent is ideal for this type of problem
    Real E = RootFinding::FindRootBrent(func, low, high, 1e-15);
    
    std::cout << "Kepler Equation (e=0.5, M=0.5):\n";
    std::cout << "  Eccentric Anomaly E = " << std::setprecision(15) << E << std::endl;
    std::cout << "  True value:          " << ROOT_KEPLER_05_05 << std::endl;
    std::cout << "  Error:               " << std::abs(E - ROOT_KEPLER_05_05) << std::endl;
}
```

---

## MML Root Finding Solver Reference

### Available Methods

| Method | Derivative? | Convergence | Bracket? | Use Case |
|--------|------------|-------------|----------|----------|
| `FindRootBisection` | No | Linear | Yes | Robust baseline |
| `FindRootNewton` | Numerical | Quadratic | Yes | Speed (smooth functions) |
| `FindRootSecant` | No | Superlinear (φ≈1.618) | No | No derivatives |
| `FindRootFalsePosition` | No | Superlinear | Yes | Conservative |
| `FindRootRidders` | No | Quadratic | Yes | Excellent general purpose |
| `FindRootBrent` | No | Superlinear (~1.839) | Yes | **Gold standard** |

### Bracketing Utilities

| Method | Purpose |
|--------|---------|
| `BracketRoot` | Expand interval until sign change found |
| `FindRootBrackets` | Find all brackets in an interval |

### Algorithm Selection Guide

```
Need guaranteed convergence?
├── Yes → Use bracket-based method
│   ├── Simple/robust? → Bisection
│   ├── General purpose? → Brent (recommended)
│   └── Quadratic without derivatives? → Ridders
│
└── No → Use non-bracket method
    ├── Have good initial guess? → Secant
    └── Have derivative? → Newton
```

---

## Test Function Categories by Difficulty

| Difficulty | Description | Examples |
|------------|-------------|----------|
| **Easy (1)** | Simple roots, well-behaved | Quadratics, Dottie, Omega |
| **Medium (2)** | Multiple roots, wider brackets | Quartic, Wilkinson, Kepler |
| **Hard (3)** | Close roots, steep functions | Close roots, Flat cubic, Atan |
| **Challenging (4)** | Very close roots, ill-conditioned | Roots separated by 0.0001 |

---

## Convergence Notes

### Newton's Method on Multiple Roots

For a root of multiplicity $m$, Newton's method converges only **linearly** instead of quadratically:

| Multiplicity | Convergence | Iterations for 10⁻¹⁰ |
|--------------|-------------|----------------------|
| m = 1 | Quadratic | ~5 |
| m = 2 | Linear | ~35 |
| m = 3 | Linear (slower) | ~50+ |

### Brent vs Ridders

Both are excellent choices:
- **Brent**: Adaptive, uses secant and inverse quadratic interpolation
- **Ridders**: Simpler, uses exponential interpolation

For most problems, they perform similarly. Brent is slightly more robust.

---

## See Also

- [Optimization_testbed.md](Optimization_testbed.md) - Optimization test cases
- [Functions_testbed.md](Functions_testbed.md) - General test functions
- [RootFinding.h](../../mml/algorithms/RootFinding.h) - MML root finding algorithms