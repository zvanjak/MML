# Integration Test Bed Documentation

## Overview

The Integration Test Bed provides a comprehensive collection of test integrals for validating and benchmarking MML's numerical integration algorithms. Located in `test_data/integration_test_bed.h` and `test_data/integration_defs.h`, it offers 15 carefully selected test cases organized by difficulty and integral type.

**Key Features:**
- 15 test integrals with known analytical results
- Three categories: smooth, oscillatory, and improper integrals
- `TestIntegral` struct with metadata (difficulty, tolerances, singularity info)
- `IntegrandWrapper` class to adapt `std::function` to `IRealFunction`
- Utility functions for result verification

## Architecture

### Core Files

| File | Purpose |
|------|---------|
| [integration_test_bed.h](../../test_data/integration_test_bed.h) | Test structures, collection generators, utilities |
| [integration_defs.h](../../test_data/integration_defs.h) | Test function implementations and exact results |

### Data Structures

```cpp
namespace MML::TestBeds {

/// Classification of integral type
enum class IntegralType {
    Regular,            ///< Normal definite integral
    Oscillatory,        ///< Highly oscillatory integrand
    SemiInfinite,       ///< One limit is ±∞
    DoublyInfinite,     ///< Both limits are ±∞
    Singular,           ///< Integrand has integrable singularity
    SingularEndpoint    ///< Singularity at endpoint
};

/// Complete test case descriptor
struct TestIntegral {
    std::string name;                            ///< Descriptive name
    std::function<Real(Real)> func;              ///< The integrand f(x)
    Real a, b;                                   ///< Integration bounds
    Real exactResult;                            ///< Known analytical result
    
    // Optional
    std::function<Real(Real)> antiderivative;    ///< F(x) where F'(x) = f(x)
    bool hasAntiderivative = false;
    
    // Characteristics
    IntegralType type = IntegralType::Regular;
    bool isImproper = false;
    Real singularityLocation = 0;
    
    // Metadata
    std::string category;                        ///< "smooth", "oscillatory", "improper"
    std::string description;
    int difficulty;                              ///< 1=easy ... 4=very hard
    Real suggestedTolerance;                     ///< Reasonable tolerance
};

}
```

### IntegrandWrapper

Adapts `std::function<Real(Real)>` to MML's `IRealFunction` interface:

```cpp
class IntegrandWrapper : public IRealFunction {
private:
    std::function<Real(Real)> _func;
public:
    IntegrandWrapper(std::function<Real(Real)> f) : _func(f) {}
    Real operator()(Real x) const override { return _func(x); }
};
```

---

## Test Categories

### Category 1: Smooth Continuous Functions

These are well-behaved integrands that converge quickly with standard quadrature rules. Ideal for basic validation.

| # | Name | Integrand f(x) | Interval | Exact Result | Difficulty |
|---|------|----------------|----------|--------------|------------|
| 1 | Polynomial | x³ - 2x + 1 | [0, 2] | 2.0 | 1 |
| 2 | Gaussian | exp(-x²/2) | [-1, 1] | 1.71124... | 1 |
| 3 | SinCos | sin(x)·cos(x) | [0, π] | 0.0 | 1 |
| 4 | Arctan | 1/(1+x²) | [0, 1] | π/4 ≈ 0.7854 | 1 |
| 5 | XExpNegX | x·exp(-x) | [0, 3] | 0.8009... | 1 |

**Function Implementations:**

```cpp
// 1. Polynomial: ∫[0,2] (x³-2x+1) dx = 2
static Real Int_Polynomial(Real x) { return x*x*x - 2*x + 1; }
static Real Int_Polynomial_antideriv(Real x) { return x*x*x*x/4 - x*x + x; }

// 2. Gaussian bell curve (no closed-form antiderivative)
static Real Int_Gaussian(Real x) { return std::exp(-x*x/2); }

// 3. Trigonometric: sin(x)cos(x) = sin(2x)/2
static Real Int_SinCos(Real x) { return std::sin(x) * std::cos(x); }
static Real Int_SinCos_antideriv(Real x) { return std::sin(x)*std::sin(x)/2; }

// 4. Arctangent derivative
static Real Int_Arctan(Real x) { return 1.0 / (1 + x*x); }
static Real Int_Arctan_antideriv(Real x) { return std::atan(x); }

// 5. Exponential decay modulated by x
static Real Int_XExpNegX(Real x) { return x * std::exp(-x); }
static Real Int_XExpNegX_antideriv(Real x) { return -std::exp(-x)*(x+1); }
```

---

### Category 2: Oscillatory Functions

Rapidly oscillating integrands that require more quadrature points for accurate results. Tests the robustness of integration routines.

| # | Name | Integrand f(x) | Interval | Exact Result | Difficulty |
|---|------|----------------|----------|--------------|------------|
| 6 | FastSine | sin(10x) | [0, π] | 0.0 | 2 |
| 7 | DampedOsc | exp(-x)·sin(5x) | [0, 2π] | ≈ 0.192 | 2 |
| 8 | FresnelCos | cos(x²) | [0, √(π/2)] | √(π/8) ≈ 0.627 | 3 |
| 9 | SinProduct | sin(x)·sin(3x) | [0, π] | 0.0 | 2 |
| 10 | Chirp | sin(x²) | [0, √(2π)] | √(π/8) ≈ 0.627 | 3 |

**Function Implementations:**

```cpp
// 6. High-frequency sine: exactly 0 over full period
static Real Int_FastSine(Real x) { return std::sin(10*x); }
static Real Int_FastSine_antideriv(Real x) { return -std::cos(10*x)/10; }

// 7. Damped oscillation (exponential × sine)
static Real Int_DampedOsc(Real x) { return std::exp(-x) * std::sin(5*x); }
static Real Int_DampedOsc_antideriv(Real x) { 
    return std::exp(-x) * (-std::sin(5*x) - 5*std::cos(5*x)) / 26; 
}

// 8. Fresnel-type: cos(x²) - frequency increases with x
static Real Int_CosSq(Real x) { return std::cos(x*x); }

// 9. Product of sines (orthogonality test)
static Real Int_SinProduct(Real x) { return std::sin(x) * std::sin(3*x); }

// 10. Chirp signal: sin(x²) - no closed-form antiderivative
static Real Int_Chirp(Real x) { return std::sin(x*x); }
```

---

### Category 3: Improper Integrals

Integrals with singularities or infinite limits. Tests special handling capabilities of integration algorithms.

| # | Name | Integrand f(x) | Interval | Exact Result | Singularity |
|---|------|----------------|----------|--------------|-------------|
| 11 | SqrtSing | 1/√x | [0, 1] | 2.0 | x = 0 |
| 12 | LogSing | -ln(x) | [0, 1] | 1.0 | x = 0 |
| 13 | ExpDecay | exp(-x) | [0, 20]* | ≈ 1.0 | ∞ (truncated) |
| 14 | GaussTail | exp(-x²) | [0, 6]* | √π/2 ≈ 0.886 | ∞ (truncated) |
| 15 | AlgDecay | 1/(1+x²) | [0, 1000]* | π/2 ≈ 1.571 | ∞ (truncated) |

*Using finite cutoff for practical computation

**Function Implementations:**

```cpp
// 11. Square root singularity at origin
static Real Int_SqrtSing(Real x) { 
    if (x <= 0) return std::numeric_limits<Real>::infinity();
    return 1.0 / std::sqrt(x); 
}
static Real Int_SqrtSing_antideriv(Real x) { return 2*std::sqrt(x); }

// 12. Logarithmic singularity at origin
static Real Int_LogSing(Real x) { 
    if (x <= 0) return std::numeric_limits<Real>::infinity();
    return -std::log(x); 
}

// 13-15. Semi-infinite integrals (use finite cutoffs)
static Real Int_ExpDecay(Real x) { return std::exp(-x); }
static Real Int_GaussTail(Real x) { return std::exp(-x*x); }
static Real Int_AlgDecay(Real x) { return 1.0 / (1 + x*x); }
```

---

## API Reference

### Test Collection Generators

```cpp
namespace MML::TestBeds {

/// Get all 5 smooth test cases
std::vector<TestIntegral> getSmoothIntegrals();

/// Get all 5 oscillatory test cases  
std::vector<TestIntegral> getOscillatoryIntegrals();

/// Get all 5 improper test cases
std::vector<TestIntegral> getImproperIntegrals();

/// Get all 15 test cases
std::vector<TestIntegral> getAllIntegrals();

/// Get tests filtered by difficulty (1-4)
std::vector<TestIntegral> getIntegralsByDifficulty(int maxDifficulty);

}
```

### Individual Test Accessors

```cpp
TestIntegral getPolynomialIntegral();     // Smooth #1
TestIntegral getGaussianIntegral();       // Smooth #2
TestIntegral getSinCosIntegral();         // Smooth #3
TestIntegral getArctanIntegral();         // Smooth #4
TestIntegral getXExpNegXIntegral();       // Smooth #5

TestIntegral getFastSineIntegral();       // Oscillatory #6
TestIntegral getDampedOscIntegral();      // Oscillatory #7
TestIntegral getFresnelCosIntegral();     // Oscillatory #8
TestIntegral getSinProductIntegral();     // Oscillatory #9
TestIntegral getChirpIntegral();          // Oscillatory #10

TestIntegral getSqrtSingularityIntegral();     // Improper #11
TestIntegral getLogSingularityIntegral();      // Improper #12
TestIntegral getExpDecayInfiniteIntegral();    // Improper #13
TestIntegral getGaussianTailIntegral();        // Improper #14
TestIntegral getAlgebraicDecayIntegral();      // Improper #15
```

### Utility Functions

```cpp
/// Create MML-compatible function wrapper from test case
IntegrandWrapper makeIntegrand(const TestIntegral& test);

/// Check if computed result matches expected within tolerance
bool checkResult(Real computed, Real expected, Real tolerance);
```

---

## Usage Examples

### Basic Integration with Test Bed

```cpp
#include "test_data/integration_test_bed.h"
#include "core/Integration.h"

using namespace MML;
using namespace MML::TestBeds;

// Get a test case and integrate it
TestIntegral poly = getPolynomialIntegral();
IntegrandWrapper f = makeIntegrand(poly);

// Using Simpson's rule
IntegrationResult result = IntegrateSimpson(f, poly.a, poly.b, 1e-10);

std::cout << "Computed: " << result.value << "\n";
std::cout << "Expected: " << poly.exactResult << "\n";
std::cout << "Error: " << result.error_estimate << "\n";
std::cout << "Iterations: " << result.iterations << "\n";
std::cout << "Converged: " << (result.converged ? "yes" : "no") << "\n";
```

### Comparing Integration Methods

```cpp
#include "test_data/integration_test_bed.h"
#include "core/Integration.h"

using namespace MML;
using namespace MML::TestBeds;

TestIntegral test = getGaussianIntegral();
IntegrandWrapper f = makeIntegrand(test);
Real tol = 1e-8;

// Trapezoidal rule
IntegrationResult trapResult = IntegrateTrap(f, test.a, test.b, tol);

// Simpson's rule  
IntegrationResult simpResult = IntegrateSimpson(f, test.a, test.b, tol);

// Romberg extrapolation
IntegrationResult rombResult = IntegrateRomberg(f, test.a, test.b, tol);

// 10-point Gauss-Legendre (no tolerance parameter)
IntegrationResult gaussResult = IntegrateGauss10(f, test.a, test.b);

std::cout << "Method        Value          Iterations  Error Est\n";
std::cout << "Trapezoidal   " << trapResult.value  << "  " 
          << trapResult.iterations << "  " << trapResult.error_estimate << "\n";
std::cout << "Simpson       " << simpResult.value  << "  "
          << simpResult.iterations << "  " << simpResult.error_estimate << "\n";
std::cout << "Romberg       " << rombResult.value  << "  "
          << rombResult.iterations << "  " << rombResult.error_estimate << "\n";
std::cout << "Gauss10       " << gaussResult.value << "  "
          << gaussResult.iterations << "  (fixed-order)\n";
```

### Running All Tests in a Category

```cpp
#include "test_data/integration_test_bed.h"
#include "core/Integration.h"

using namespace MML;
using namespace MML::TestBeds;

auto tests = getSmoothIntegrals();

for (const auto& test : tests) {
    IntegrandWrapper f = makeIntegrand(test);
    IntegrationResult result = IntegrateRomberg(f, test.a, test.b, test.suggestedTolerance);
    
    bool passed = checkResult(result.value, test.exactResult, test.suggestedTolerance);
    
    std::cout << test.name << ": ";
    std::cout << (passed ? "PASS" : "FAIL");
    std::cout << " (computed=" << result.value 
              << ", expected=" << test.exactResult << ")\n";
}
```

### Handling Oscillatory Integrals

```cpp
#include "test_data/integration_test_bed.h"
#include "core/Integration.h"

using namespace MML;
using namespace MML::TestBeds;

// Oscillatory integrals often need tighter tolerance or Romberg
auto oscillatory = getOscillatoryIntegrals();

for (const auto& test : oscillatory) {
    IntegrandWrapper f = makeIntegrand(test);
    
    // Romberg works better for oscillatory functions
    IntegrationResult result = IntegrateRomberg(f, test.a, test.b, 1e-8);
    
    if (!result.converged) {
        std::cerr << "Warning: " << test.name << " did not converge!\n";
    }
    
    Real relError = std::abs(result.value - test.exactResult);
    if (test.exactResult != 0) {
        relError /= std::abs(test.exactResult);
    }
    
    std::cout << test.name << ": rel_error = " << relError << "\n";
}
```

### Using Gaussian Quadrature for Improper Integrals

```cpp
#include "test_data/integration_test_bed.h"
#include "core/Integration.h"
#include "core/Integration/GaussianQuadrature.h"

using namespace MML;
using namespace MML::TestBeds;

// For semi-infinite integrals like ∫[0,∞) f(x) dx,
// Gauss-Laguerre quadrature is ideal

// Example: ∫[0,∞) x·exp(-x) dx = 1
// Gauss-Laguerre integrates f(x)·e^(-x), so we supply f(x) = x
RealFunctionFromStdFunc f([](Real x) { return x; });
IntegrationResult result = IntegrateGaussLaguerre(f, 20, 0.0);
std::cout << "∫x·e^(-x) = " << result.value << " (exact: 1)\n";

// For doubly-infinite ∫(-∞,∞) f(x)·exp(-x²) dx, use Gauss-Hermite
// Example: ∫ e^(-x²) dx = √π
// Gauss-Hermite has weight e^(-x²) built in, so f(x) = 1
RealFunctionFromStdFunc g([](Real x) { return 1.0; });
IntegrationResult result2 = IntegrateGaussHermite(g, 20);
std::cout << "∫e^(-x²) = " << result2.value << " (exact: " << std::sqrt(Constants::PI) << ")\n";
```

### Filtering by Difficulty

```cpp
#include "test_data/integration_test_bed.h"

using namespace MML::TestBeds;

// Get only "easy" tests (difficulty 1)
auto easyTests = getIntegralsByDifficulty(1);
std::cout << "Easy tests: " << easyTests.size() << "\n";  // 5

// Get easy and medium tests (difficulty <= 2)
auto normalTests = getIntegralsByDifficulty(2);
std::cout << "Normal tests: " << normalTests.size() << "\n";  // 10

// Get all tests including hard ones
auto allTests = getIntegralsByDifficulty(4);
std::cout << "All tests: " << allTests.size() << "\n";  // 15
```

### Using the Unified Integrate Function

```cpp
#include "test_data/integration_test_bed.h"
#include "core/Integration.h"

using namespace MML;
using namespace MML::TestBeds;

TestIntegral test = getArctanIntegral();
IntegrandWrapper f = makeIntegrand(test);

// Explicitly specify integration method
IntegrationResult r1 = Integrate(f, test.a, test.b, TRAP, 1e-8);
IntegrationResult r2 = Integrate(f, test.a, test.b, SIMPSON, 1e-8);
IntegrationResult r3 = Integrate(f, test.a, test.b, ROMBERG, 1e-10);
IntegrationResult r4 = Integrate(f, test.a, test.b, GAUSS10);  // No tolerance

std::cout << "TRAP:    " << r1.value << " (iter=" << r1.iterations << ")\n";
std::cout << "SIMPSON: " << r2.value << " (iter=" << r2.iterations << ")\n";
std::cout << "ROMBERG: " << r3.value << " (iter=" << r3.iterations << ")\n";
std::cout << "GAUSS10: " << r4.value << " (fixed 10 pts)\n";
```

---

## MML Integration API Summary

### Core Functions (all return `IntegrationResult`)

| Function | Description | Parameters |
|----------|-------------|------------|
| `IntegrateTrap(f, a, b, eps)` | Adaptive trapezoidal | eps = tolerance (default 1e-6) |
| `IntegrateSimpson(f, a, b, eps)` | Adaptive Simpson's rule | eps = tolerance (default 1e-6) |
| `IntegrateRomberg(f, a, b, eps)` | Romberg extrapolation | eps = tolerance (default 1e-8) |
| `IntegrateGauss10(f, a, b)` | 10-point Gauss-Legendre | No tolerance (fixed-order) |
| `Integrate(f, a, b, method, eps)` | Unified interface | method = TRAP/SIMPSON/ROMBERG/GAUSS10 |

### Gaussian Quadrature (all return `IntegrationResult`)

| Function | Domain | Weight Function |
|----------|--------|-----------------|
| `IntegrateGaussLegendre(f, a, b, n)` | [a, b] | 1 (no weight) |
| `IntegrateGaussLaguerre(f, n, α)` | [0, ∞) | x^α · e^(-x) |
| `IntegrateGaussHermite(f, n)` | (-∞, ∞) | e^(-x²) |
| `IntegrateGaussJacobi(f, α, β, n)` | [-1, 1] | (1-x)^α · (1+x)^β |
| `IntegrateGaussChebyshev1(f, n)` | [-1, 1] | 1/√(1-x²) |
| `IntegrateGaussChebyshev2(f, n)` | [-1, 1] | √(1-x²) |

### IntegrationResult Structure

```cpp
struct IntegrationResult {
    Real value;           ///< Computed integral value
    Real error_estimate;  ///< Estimated absolute error
    int iterations;       ///< Number of iterations/refinements
    bool converged;       ///< True if tolerance was achieved
    
    operator Real() const;  ///< Implicit conversion (deprecated)
};
```

---

## Test Design Philosophy

1. **Known Analytical Results**: Every test has an exact or high-precision reference value
2. **Graded Difficulty**: Tests range from trivial (polynomials) to challenging (Fresnel integrals)
3. **Category Coverage**: Tests smooth, oscillatory, and singular integrands
4. **Suggested Tolerances**: Each test includes a reasonable tolerance target
5. **Rich Metadata**: Descriptions explain what each test validates

---

## See Also

- [Integration.md](../core/Integration.md) - Core integration algorithms documentation
- [GaussianQuadrature.h](../../mml/core/Integration/GaussianQuadrature.h) - Gaussian quadrature implementation
- [real_functions_test_bed.h](../../test_data/real_functions_test_bed.h) - Related function test bed