# Function Analyzer

**Location**: `mml/algorithms/FunctionsAnalyzer.h`  
**Dependencies**: RootFinding, Derivation, Integration, Statistics

---

## Overview

The Function Analyzer framework provides comprehensive tools for analyzing real-valued functions: finding roots, extrema, inflection points, detecting discontinuities, and comparing functions. It combines numerical techniques from calculus, root finding, and statistics to provide deep insights into function behavior.

### Key Capabilities

- **Critical Point Analysis**: Locate and classify local minima, maxima, and saddle points
- **Discontinuity Detection**: Find and classify jump, removable, infinite, and oscillatory discontinuities
- **Continuity Verification**: Test function continuity at points and over intervals
- **Root Finding**: Locate zeros of functions in intervals
- **Inflection Points**: Find where concavity changes (second derivative sign changes)
- **Monotonicity Testing**: Determine if functions are strictly increasing/decreasing
- **Function Comparison**: Compute difference metrics (absolute, relative, integrated)
- **Limit Computation**: Calculate left and right limits numerically

### Applications

- **Engineering**: Signal analysis, stability studies, control systems
- **Physics**: Phase transitions, equilibrium points, potential wells
- **Economics**: Optimization, elasticity analysis, marginal analysis
- **Data Science**: Model validation, function fitting assessment
- **Mathematics**: Calculus education, function classification, behavior study

---

## Quick Reference

### Analysis Classes

| Class | Purpose | Key Methods |
|-------|---------|-------------|
| `RealFunctionAnalyzer` | Single function analysis | `GetRoots`, `GetLocalOptimums`, `FindDiscontinuities` |
| `RealFunctionComparer` | Compare two functions | `getAbsDiffMax`, `getIntegratedAbsDiff` |
| `ScalarFieldAnalyzer` | 3D scalar field analysis | `IsConservative` (gradient field check) |
| `VectorFieldAnalyzer` | 3D vector field analysis | `IsSolenoidal` (divergence-free check) |

### Discontinuity Types

| Type | Description | Example |
|------|-------------|---------|
| `JUMP` | Left/right limits differ | Step function, `floor(x)` |
| `REMOVABLE` | Limits match but f(x) differs/undefined | `sin(x)/x` at x=0 |
| `INFINITE` | At least one limit is infinite | `1/x` at x=0 |
| `OSCILLATORY` | Limit doesn't exist due to oscillation | `sin(1/x)` at x=0 |

### Critical Point Types

| Type | Condition | f'(x) | f''(x) |
|------|-----------|-------|--------|
| `LOCAL_MINIMUM` | Valley | ≈ 0 | > 0 |
| `LOCAL_MAXIMUM` | Peak | ≈ 0 | < 0 |
| `SADDLE_POINT` | Inflection with f'=0 | ≈ 0 | ≈ 0 (changes sign) |

---

## RealFunctionAnalyzer - Core Analysis

### Constructor

```cpp
class RealFunctionAnalyzer
{
public:
    RealFunctionAnalyzer(const IRealFunction& f);
    RealFunctionAnalyzer(const IRealFunction& f, std::string funcName);
};
```

**Parameters:**
- `f` - Function implementing `IRealFunction` interface
- `funcName` - Optional descriptive name for output formatting

**Example 1: Basic Setup**

```cpp
#include "algorithms/FunctionsAnalyzer.h"

// Polynomial: f(x) = x³ - 6x² + 9x + 1
auto polynomial = [](Real x) { return x*x*x - 6*x*x + 9*x + 1; };
RealFunctionFromStdFunc func(polynomial);

RealFunctionAnalyzer analyzer(func, "x³ - 6x² + 9x + 1");
```

---

## Root Finding

### GetRoots - Find All Zeros in Interval

```cpp
std::vector<Real> GetRoots(Real x1, Real x2, Real eps = 1e-6);
```

**Algorithm**:
1. Create derivative wrapper f'(x)
2. Find roots of derivative (critical points)
3. Use root bracketing to find all sign changes
4. Refine with bisection to tolerance `eps`

**Example 2: Finding Roots**

```cpp
// Function: f(x) = x³ - 2x - 5
auto cubic = [](Real x) { return x*x*x - 2*x - 5; };
RealFunctionFromStdFunc func(cubic);
RealFunctionAnalyzer analyzer(func);

// Find all roots in [-10, 10]
std::vector<Real> roots = analyzer.GetRoots(-10.0, 10.0, 1e-8);

for (Real root : roots)
{
    std::cout << "Root found at x = " << root 
              << ", f(x) = " << func(root) << std::endl;
}
// Output: Root found at x ≈ 2.0946, f(x) ≈ 0
```

**Example 3: Periodic Function Roots**

```cpp
// Sine function with phase shift
auto sineFunc = [](Real x) { return std::sin(x - 0.5); };
RealFunctionFromStdFunc func(sineFunc);
RealFunctionAnalyzer analyzer(func);

// Find roots in [0, 4π]
std::vector<Real> roots = analyzer.GetRoots(0.0, 4.0 * Constants::PI);

// Results: roots at 0.5, π+0.5, 2π+0.5, 3π+0.5
std::cout << "Found " << roots.size() << " roots" << std::endl;
for (Real root : roots)
    std::cout << "  x = " << root << std::endl;
```

### calcRootsPeriod - Average Period Between Roots

```cpp
Real calcRootsPeriod(Real t1, Real t2, int numPoints);
```

**Use Case**: Detect periodicity in oscillating functions.

**Example 4: Period Detection**

```cpp
// Damped oscillator: sin(5x) * exp(-0.1x)
auto dampedOsc = [](Real x) { 
    return std::sin(5.0 * x) * std::exp(-0.1 * x); 
};
RealFunctionFromStdFunc func(dampedOsc);
RealFunctionAnalyzer analyzer(func);

// Calculate average period between zeros
Real period = analyzer.calcRootsPeriod(0.0, 20.0, 1000);
std::cout << "Average period: " << period << std::endl;
// Expected: ≈ 2π/5 ≈ 1.257 (frequency = 5)
```

---

## Critical Points - Extrema Analysis

### GetLocalOptimums - Find All Local Min/Max

```cpp
Vector<Real> GetLocalOptimums(Real x1, Real x2, Real eps = 1e-6);
```

**Algorithm**:
1. Find roots of f'(x) = 0
2. For each critical point, verify f''(x) ≠ 0
3. Classify using second derivative test

**Example 5: Find All Extrema**

```cpp
// Polynomial with multiple extrema: f(x) = x⁴ - 4x³ + 4x²
auto poly = [](Real x) { 
    return x*x*x*x - 4*x*x*x + 4*x*x; 
};
RealFunctionFromStdFunc func(poly);
RealFunctionAnalyzer analyzer(func);

Vector<Real> extrema = analyzer.GetLocalOptimums(-2.0, 4.0);

std::cout << "Found " << extrema.size() << " local optimums:" << std::endl;
for (Real x : extrema)
{
    Real value = func(x);
    Real secDeriv = Derivation::NSecDer4(func, x);
    
    std::cout << "  x = " << x << ", f(x) = " << value;
    if (secDeriv > 0)
        std::cout << " (local minimum)" << std::endl;
    else
        std::cout << " (local maximum)" << std::endl;
}
```

### GetLocalOptimumsClassified - Extrema with Type

```cpp
std::vector<CriticalPoint> GetLocalOptimumsClassified(Real x1, Real x2, Real eps = 1e-6);
```

**Returns**: Vector of `CriticalPoint` structs with classification.

**CriticalPoint Structure**:
```cpp
struct CriticalPoint
{
    Real x;                       // Location
    CriticalPointType type;       // LOCAL_MINIMUM, LOCAL_MAXIMUM, SADDLE_POINT
    Real value;                   // f(x) at critical point
};
```

**Example 6: Classified Critical Points**

```cpp
// Function with multiple features: sin(x) + 0.5*sin(3x)
auto composite = [](Real x) { 
    return std::sin(x) + 0.5 * std::sin(3.0 * x); 
};
RealFunctionFromStdFunc func(composite);
RealFunctionAnalyzer analyzer(func, "sin(x) + 0.5*sin(3x)");

auto criticalPoints = analyzer.GetLocalOptimumsClassified(0.0, 2*Constants::PI);

for (const auto& cp : criticalPoints)
{
    std::cout << "Critical point at x = " << cp.x << std::endl;
    std::cout << "  Value: f(x) = " << cp.value << std::endl;
    std::cout << "  Type: ";
    
    switch (cp.type)
    {
        case CriticalPointType::LOCAL_MINIMUM:
            std::cout << "Local MINIMUM" << std::endl;
            break;
        case CriticalPointType::LOCAL_MAXIMUM:
            std::cout << "Local MAXIMUM" << std::endl;
            break;
        case CriticalPointType::SADDLE_POINT:
            std::cout << "SADDLE POINT" << std::endl;
            break;
    }
}
```

### MinInNPoints / MaxInNPoints - Sampling Methods

```cpp
Real MinInNPoints(Real x1, Real x2, int numPoints);
Real MaxInNPoints(Real x1, Real x2, int numPoints);
```

**Use Case**: Quick min/max estimation via uniform sampling (not guaranteed global).

**Example 7: Quick Min/Max**

```cpp
RealFunctionAnalyzer analyzer(func);

// Sample 1000 points to find approximate min/max
Real minVal = analyzer.MinInNPoints(0.0, 10.0, 1000);
Real maxVal = analyzer.MaxInNPoints(0.0, 10.0, 1000);

std::cout << "Approximate minimum: " << minVal << std::endl;
std::cout << "Approximate maximum: " << maxVal << std::endl;
```

---

## Inflection Points

### GetInflectionPoints - Find Concavity Changes

```cpp
Vector<Real> GetInflectionPoints(Real x1, Real x2, Real eps = 1e-6);
```

**Algorithm**:
1. Find roots of f''(x) = 0
2. Verify second derivative sign changes (concavity change)
3. Return all confirmed inflection points

**Example 8: Inflection Points**

```cpp
// Cubic with inflection: f(x) = x³ - 3x² + 2x
auto cubic = [](Real x) { return x*x*x - 3*x*x + 2*x; };
RealFunctionFromStdFunc func(cubic);
RealFunctionAnalyzer analyzer(func);

Vector<Real> inflections = analyzer.GetInflectionPoints(-2.0, 4.0);

for (Real x : inflections)
{
    std::cout << "Inflection point at x = " << x << std::endl;
    std::cout << "  f(x) = " << func(x) << std::endl;
    
    // Verify concavity change
    Real leftConcav = Derivation::NSecDer4(func, x - 0.01);
    Real rightConcav = Derivation::NSecDer4(func, x + 0.01);
    
    std::cout << "  f''(x-ε) = " << leftConcav << std::endl;
    std::cout << "  f''(x+ε) = " << rightConcav << std::endl;
    std::cout << "  Concavity changes: " << (leftConcav * rightConcav < 0 ? "YES" : "NO") << std::endl;
}
// Expected: inflection at x = 1.0
```

---

## Continuity Analysis

### Point-Level Tests

```cpp
bool isDefinedAtPoint(Real x) const;
bool isContinuousAtPoint(Real x, Real eps = 1e-6);
bool isDerivativeDefinedAtPoint(Real x, Real eps = 1e-6);
bool isInflectionPoint(Real x, Real eps = 1e-6);
bool isLocalOptimum(Real x, Real eps = 1e-6);
```

**Example 9: Point Analysis**

```cpp
// Piecewise function with discontinuity at x=0
auto piecewise = [](Real x) {
    return (x >= 0) ? x*x : -x;
};
RealFunctionFromStdFunc func(piecewise);
RealFunctionAnalyzer analyzer(func);

// Test at discontinuity
Real x = 0.0;
std::cout << "Analysis at x = " << x << ":" << std::endl;
std::cout << "  Defined: " << (analyzer.isDefinedAtPoint(x) ? "yes" : "no") << std::endl;
std::cout << "  Continuous: " << (analyzer.isContinuousAtPoint(x, 1e-6) ? "yes" : "no") << std::endl;
std::cout << "  Derivative exists: " << (analyzer.isDerivativeDefinedAtPoint(x, 1e-6) ? "yes" : "no") << std::endl;

// Output:
// Defined: yes
// Continuous: yes
// Derivative exists: no (left deriv = -1, right deriv = 0)
```

### Interval-Level Tests

```cpp
bool isContinuous(Real x1, Real x2, int numPoints, Real eps = 1e-6, 
                  std::vector<Real>* discontinuities = nullptr);

bool isMonotonic(Real x1, Real x2, int numPoints);
```

**Example 10: Interval Continuity**

```cpp
// Step function
auto stepFunc = [](Real x) { 
    return std::floor(x); 
};
RealFunctionFromStdFunc func(stepFunc);
RealFunctionAnalyzer analyzer(func);

// Check continuity with discontinuity detection
std::vector<Real> discontinuities;
bool continuous = analyzer.isContinuous(0.0, 5.0, 500, 1e-6, &discontinuities);

std::cout << "Function is " << (continuous ? "continuous" : "discontinuous") << std::endl;

if (!discontinuities.empty())
{
    std::cout << "Discontinuities found at:" << std::endl;
    for (Real x : discontinuities)
        std::cout << "  x ≈ " << x << std::endl;
    // Expected: discontinuities near 1, 2, 3, 4
}
```

---

## Discontinuity Detection and Classification

### Limit Computation

```cpp
Real ComputeLeftLimit(Real x, Real eps = 1e-6, int maxSteps = 20) const;
Real ComputeRightLimit(Real x, Real eps = 1e-6, int maxSteps = 20) const;
```

**Algorithm**: Approach point with decreasing step sizes, check convergence.

**Example 11: Computing Limits**

```cpp
// Function with jump: f(x) = (x < 0) ? -1 : 1
auto signFunc = [](Real x) { return (x < 0) ? -1.0 : 1.0; };
RealFunctionFromStdFunc func(signFunc);
RealFunctionAnalyzer analyzer(func);

Real x = 0.0;
Real leftLimit = analyzer.ComputeLeftLimit(x);
Real rightLimit = analyzer.ComputeRightLimit(x);

std::cout << "At x = " << x << ":" << std::endl;
std::cout << "  lim(x→0⁻) f(x) = " << leftLimit << std::endl;    // -1
std::cout << "  lim(x→0⁺) f(x) = " << rightLimit << std::endl;   // +1
std::cout << "  f(0) = " << func(0) << std::endl;                 // +1
std::cout << "  Jump size: " << std::abs(rightLimit - leftLimit) << std::endl;  // 2
```

### ClassifyDiscontinuity - Detailed Classification

```cpp
DiscontinuityPoint ClassifyDiscontinuity(Real x, Real eps = 1e-6) const;
```

**DiscontinuityPoint Structure**:
```cpp
struct DiscontinuityPoint
{
    Real x;                       // Location
    DiscontinuityType type;       // JUMP, REMOVABLE, INFINITE, OSCILLATORY, UNKNOWN
    Real leftLimit;               // lim(x→a⁻) f(x)
    Real rightLimit;              // lim(x→a⁺) f(x)
    Real valueAtPoint;            // f(a)
    Real jumpSize;                // |rightLimit - leftLimit|
};
```

**Example 12: Classify Jump Discontinuity**

```cpp
// Heaviside step function
auto heaviside = [](Real x) { return (x >= 0) ? 1.0 : 0.0; };
RealFunctionFromStdFunc func(heaviside);
RealFunctionAnalyzer analyzer(func);

DiscontinuityPoint disc = analyzer.ClassifyDiscontinuity(0.0);

std::cout << "Discontinuity at x = " << disc.x << std::endl;
std::cout << "Type: ";
switch (disc.type)
{
    case DiscontinuityType::JUMP:
        std::cout << "JUMP (size: " << disc.jumpSize << ")" << std::endl;
        break;
    case DiscontinuityType::REMOVABLE:
        std::cout << "REMOVABLE" << std::endl;
        break;
    case DiscontinuityType::INFINITE:
        std::cout << "INFINITE" << std::endl;
        break;
    default:
        std::cout << "OTHER" << std::endl;
}

std::cout << "Left limit:  " << disc.leftLimit << std::endl;
std::cout << "Right limit: " << disc.rightLimit << std::endl;
std::cout << "f(0):        " << disc.valueAtPoint << std::endl;
```

**Example 13: Removable Discontinuity**

```cpp
// Classic removable: f(x) = sin(x)/x, undefined at x=0
auto sinc = [](Real x) { 
    return (x == 0.0) ? std::numeric_limits<Real>::quiet_NaN() : std::sin(x) / x; 
};
RealFunctionFromStdFunc func(sinc);
RealFunctionAnalyzer analyzer(func);

DiscontinuityPoint disc = analyzer.ClassifyDiscontinuity(0.0);

std::cout << "Type: " << (disc.type == DiscontinuityType::REMOVABLE ? "REMOVABLE" : "OTHER") << std::endl;
std::cout << "Left limit:  " << disc.leftLimit << std::endl;   // ≈ 1.0
std::cout << "Right limit: " << disc.rightLimit << std::endl;  // ≈ 1.0
std::cout << "f(0):        " << (std::isnan(disc.valueAtPoint) ? "undefined" : std::to_string(disc.valueAtPoint)) << std::endl;
```

**Example 14: Infinite Discontinuity**

```cpp
// Vertical asymptote: f(x) = 1/x
auto reciprocal = [](Real x) { return 1.0 / x; };
RealFunctionFromStdFunc func(reciprocal);
RealFunctionAnalyzer analyzer(func);

DiscontinuityPoint disc = analyzer.ClassifyDiscontinuity(0.0);

std::cout << "Type: " << (disc.type == DiscontinuityType::INFINITE ? "INFINITE" : "OTHER") << std::endl;
std::cout << "Left limit:  " << disc.leftLimit << std::endl;   // -∞
std::cout << "Right limit: " << disc.rightLimit << std::endl;  // +∞
```

### FindDiscontinuities - Comprehensive Search

```cpp
std::vector<DiscontinuityPoint> FindDiscontinuities(Real x1, Real x2, 
                                                     int numPoints = 1000, 
                                                     Real eps = 1e-6);
```

**Multi-Stage Algorithm**:
1. Sample at `numPoints` locations, test continuity
2. Check for rapid value changes between samples
3. Refine search in suspicious intervals
4. Classify each detected discontinuity
5. Return sorted list

**Example 15: Find All Discontinuities**

```cpp
// Tan function with vertical asymptotes
auto tanFunc = [](Real x) { return std::tan(x); };
RealFunctionFromStdFunc func(tanFunc);
RealFunctionAnalyzer analyzer(func);

// Find discontinuities in [-π, π]
auto discontinuities = analyzer.FindDiscontinuities(-Constants::PI, Constants::PI, 2000);

std::cout << "Found " << discontinuities.size() << " discontinuities:" << std::endl;
for (const auto& disc : discontinuities)
{
    std::cout << "\nAt x ≈ " << disc.x << ":" << std::endl;
    std::cout << "  Type: ";
    switch (disc.type)
    {
        case DiscontinuityType::INFINITE:
            std::cout << "INFINITE (vertical asymptote)" << std::endl;
            break;
        default:
            std::cout << "OTHER" << std::endl;
    }
}
// Expected: discontinuities near -π/2, π/2
```

### PrintContinuityAnalysis - Formatted Report

```cpp
void PrintContinuityAnalysis(Real x1, Real x2, int numPoints = 1000, Real eps = 1e-6);
```

**Example 16: Full Continuity Report**

```cpp
// Piecewise function with multiple discontinuities
auto complex = [](Real x) {
    if (x < -1) return x + 2;
    if (x < 0) return 1.0 / x;  // Infinite discontinuity at x=0
    if (x < 1) return 0.0;      // Jump at x=0
    return x - 1;               // Jump at x=1
};
RealFunctionFromStdFunc func(complex);
RealFunctionAnalyzer analyzer(func);

// Generate comprehensive analysis report
analyzer.PrintContinuityAnalysis(-2.0, 2.0, 1000);

/* Output:
=== CONTINUITY ANALYSIS ===
Interval: [-2.000000, 2.000000]
Sample points: 1000
Tolerance: 0.000001

✗ Function has 2 discontinuity point(s):

Discontinuity #1 at x = 0.000000
  Type: INFINITE
  Left limit:  -∞
  Right limit: 0.000000
  f(0.000000): 0.000000

Discontinuity #2 at x = 1.000000
  Type: JUMP (jump size: 1.000000)
  Left limit:  0.000000
  Right limit: 0.000000
  f(1.000000): 0.000000
*/
```

---

## Function Comparison - RealFunctionComparer

### Constructor

```cpp
class RealFunctionComparer
{
public:
    RealFunctionComparer(IRealFunction& f1, IRealFunction& f2);
};
```

### Absolute Difference Metrics

```cpp
Real getAbsDiffSum(Real a, Real b, int numPoints);
Real getAbsDiffAvg(Real a, Real b, int numPoints);
Real getAbsDiffMax(Real a, Real b, int numPoints);
```

**Formula**: |f₁(x) - f₂(x)|

**Example 17: Compare Approximations**

```cpp
// Original function
auto original = [](Real x) { return std::exp(x); };
RealFunctionFromStdFunc f1(original);

// Taylor approximation: exp(x) ≈ 1 + x + x²/2 + x³/6
auto taylor3 = [](Real x) { 
    return 1.0 + x + 0.5*x*x + x*x*x/6.0; 
};
RealFunctionFromStdFunc f2(taylor3);

RealFunctionComparer comparer(f1, f2);

// Compare in [-1, 1]
Real maxDiff = comparer.getAbsDiffMax(-1.0, 1.0, 1000);
Real avgDiff = comparer.getAbsDiffAvg(-1.0, 1.0, 1000);

std::cout << "Maximum absolute difference: " << maxDiff << std::endl;
std::cout << "Average absolute difference: " << avgDiff << std::endl;
```

### Relative Difference Metrics

```cpp
Real getRelDiffSum(Real a, Real b, int numPoints);
Real getRelDiffAvg(Real a, Real b, int numPoints);
Real getRelDiffMax(Real a, Real b, int numPoints);
```

**Formula**: |f₁(x) - f₂(x)| / |f₁(x)|

**Example 18: Relative Error Analysis**

```cpp
// True sine vs 3rd order approximation
auto trueSine = [](Real x) { return std::sin(x); };
auto sineApprox = [](Real x) { return x - x*x*x/6.0; };

RealFunctionFromStdFunc f1(trueSine);
RealFunctionFromStdFunc f2(sineApprox);

RealFunctionComparer comparer(f1, f2);

// Relative error in [-π/4, π/4]
Real maxRelDiff = comparer.getRelDiffMax(-Constants::PI/4, Constants::PI/4, 500);
Real avgRelDiff = comparer.getRelDiffAvg(-Constants::PI/4, Constants::PI/4, 500);

std::cout << "Maximum relative error: " << maxRelDiff * 100 << "%" << std::endl;
std::cout << "Average relative error: " << avgRelDiff * 100 << "%" << std::endl;
```

### Integrated Difference Metrics

```cpp
Real getIntegratedDiff(Real a, Real b, IntegrationMethod method = IntegrationMethod::TRAP);
Real getIntegratedAbsDiff(Real a, Real b, IntegrationMethod method = IntegrationMethod::TRAP);
Real getIntegratedSqrDiff(Real a, Real b, IntegrationMethod method = IntegrationMethod::TRAP);
```

**Formulas**:
- Integrated diff: ∫[f₁(x) - f₂(x)]dx
- Integrated abs diff: ∫|f₁(x) - f₂(x)|dx
- Integrated squared diff: ∫[f₁(x) - f₂(x)]²dx (L² norm)

**Example 19: L² Norm Distance**

```cpp
// Compare two bell curves
auto gaussian1 = [](Real x) { return std::exp(-x*x); };
auto gaussian2 = [](Real x) { return std::exp(-(x-0.5)*(x-0.5)); };  // Shifted

RealFunctionFromStdFunc f1(gaussian1);
RealFunctionFromStdFunc f2(gaussian2);

RealFunctionComparer comparer(f1, f2);

// L² distance between functions
Real sqrDiff = comparer.getIntegratedSqrDiff(-5.0, 5.0, IntegrationMethod::SIMPSON);
Real l2Distance = std::sqrt(sqrDiff);

std::cout << "L² distance: " << l2Distance << std::endl;

// Also compute absolute integrated difference
Real absDiff = comparer.getIntegratedAbsDiff(-5.0, 5.0, IntegrationMethod::SIMPSON);
std::cout << "Integrated absolute difference: " << absDiff << std::endl;
```

---

## Print Methods - Formatted Output

### PrintPointAnalysis

```cpp
void PrintPointAnalysis(Real x, Real eps = 1e-6);
```

**Example 20: Point Report**

```cpp
auto func = [](Real x) { return std::abs(x); };
RealFunctionFromStdFunc f(func);
RealFunctionAnalyzer analyzer(f, "|x|");

analyzer.PrintPointAnalysis(0.0);

/* Output:
Function analysis at point: 0:
  Defined at point:      yes
  Continuous at point:   yes
  Derivative defined:    no
  Is inflection point:   no
*/
```

### PrintIntervalAnalysis

```cpp
void PrintIntervalAnalysis(Real x1, Real x2, int numPoints, Real eps = 1e-6);
```

**Example 21: Interval Report**

```cpp
auto func = [](Real x) { return 1.0 / (1.0 + x*x); };  // Lorentzian
RealFunctionFromStdFunc f(func);
RealFunctionAnalyzer analyzer(f, "1/(1+x²)");

analyzer.PrintIntervalAnalysis(-5.0, 5.0, 1000);

/* Output:
f(x) = 1/(1+x²) - Function analysis in interval [-5.00, 5.00] with 1000 points:
  Defined    : yes
  Continuous : yes
*/
```

---

## Advanced Features

### ScalarFieldAnalyzer - 3D Field Analysis

```cpp
class ScalarFieldAnalyzer
{
public:
    ScalarFieldAnalyzer(IScalarFunction<3>& f);
    bool IsConservative();
};
```

**Purpose**: Check if gradient field is conservative (curl = 0).

**Note**: Currently under development - checks curl of gradient field.

### VectorFieldAnalyzer - 3D Vector Analysis

```cpp
class VectorFieldAnalyzer
{
public:
    VectorFieldAnalyzer(IVectorFunction<3>& f);
    bool IsSolenoidal();
};
```

**Purpose**: Check if vector field is divergence-free (∇·F = 0).

**Note**: Currently under development.

---

## Integration with Other Modules

### With Root Finding

```cpp
#include "algorithms/RootFinding.h"
#include "algorithms/FunctionsAnalyzer.h"

// Function analyzer uses RootFinding internally
RealFunctionAnalyzer analyzer(func);
std::vector<Real> roots = analyzer.GetRoots(a, b);  // Uses RootFinding::FindRootBrackets + Bisection
```

### With Derivation

```cpp
#include "core/Derivation.h"

// All derivative-based methods use Derivation module
Real firstDeriv = Derivation::NDer4(func, x);
Real secondDeriv = Derivation::NSecDer4(func, x);

// Used internally for critical points and inflection detection
```

### With Integration

```cpp
#include "core/Integration.h"

// Function comparer uses integration for L² norms
RealFunctionComparer comparer(f1, f2);
Real sqrDiff = comparer.getIntegratedSqrDiff(a, b, IntegrationMethod::SIMPSON);
```

### With Statistics

```cpp
#include "algorithms/Statistics.h"

// Period calculation uses statistics
Real period = analyzer.calcRootsPeriod(t1, t2, numPoints);  // Uses Statistics::Avg()
```

---

## Best Practices

### Choosing Tolerance (eps)

**Continuity tests**:
- Smooth functions: `eps = 1e-6` (default)
- Noisy/irregular functions: `eps = 1e-4` or higher
- High-precision needs: `eps = 1e-10`

**Root finding**:
- Engineering accuracy: `eps = 1e-6`
- High precision: `eps = 1e-10`
- Quick estimates: `eps = 1e-4`

### Sampling Density (numPoints)

**Discontinuity detection**:
- Smooth functions: 100-500 points
- Functions with sharp features: 1000-5000 points
- Pathological cases: 10000+ points

**Root finding**:
- Well-behaved: 100-200 points for bracketing
- Oscillatory: 500-2000 points

### Performance Considerations

**Expensive operations** (in order):
1. `FindDiscontinuities` - multi-stage refinement
2. `GetLocalOptimumsClassified` - multiple derivative evaluations
3. `isContinuous` - adaptive sampling
4. `GetRoots` - root bracketing + refinement
5. `GetInflectionPoints` - second derivative roots

**Optimization tips**:
- Cache function values when possible
- Use lower `numPoints` for initial exploration
- Refine locally with higher density as needed
- Reuse `RealFunctionAnalyzer` for same function

### Numerical Stability

**Avoid**:
- Testing discontinuities at exact singularities (use nearby points)
- Very large intervals with few sample points
- Functions with extreme dynamic range

**Use**:
- Adaptive eps based on function scale
- Multiple passes with increasing resolution
- Verification of results with different methods

---

## Common Patterns

### Pattern 1: Complete Function Analysis

```cpp
void AnalyzeFunction(const IRealFunction& func, Real a, Real b)
{
    RealFunctionAnalyzer analyzer(func, "f(x)");
    
    // Find all interesting features
    auto roots = analyzer.GetRoots(a, b);
    auto extrema = analyzer.GetLocalOptimumsClassified(a, b);
    auto inflections = analyzer.GetInflectionPoints(a, b);
    auto discontinuities = analyzer.FindDiscontinuities(a, b);
    
    // Print comprehensive report
    std::cout << "=== FUNCTION ANALYSIS ===" << std::endl;
    std::cout << "Roots: " << roots.size() << std::endl;
    std::cout << "Local extrema: " << extrema.size() << std::endl;
    std::cout << "Inflection points: " << inflections.size() << std::endl;
    std::cout << "Discontinuities: " << discontinuities.size() << std::endl;
    
    analyzer.PrintContinuityAnalysis(a, b);
}
```

### Pattern 2: Approximation Quality Assessment

```cpp
Real AssessApproximation(IRealFunction& exact, IRealFunction& approx, Real a, Real b)
{
    RealFunctionComparer comparer(exact, approx);
    
    // Multiple error metrics
    Real maxAbsError = comparer.getAbsDiffMax(a, b, 1000);
    Real avgAbsError = comparer.getAbsDiffAvg(a, b, 1000);
    Real maxRelError = comparer.getRelDiffMax(a, b, 1000);
    Real l2Error = std::sqrt(comparer.getIntegratedSqrDiff(a, b, IntegrationMethod::SIMPSON));
    
    std::cout << "Approximation quality:" << std::endl;
    std::cout << "  Max absolute error: " << maxAbsError << std::endl;
    std::cout << "  Avg absolute error: " << avgAbsError << std::endl;
    std::cout << "  Max relative error: " << maxRelError * 100 << "%" << std::endl;
    std::cout << "  L² error: " << l2Error << std::endl;
    
    return l2Error;
}
```

### Pattern 3: Adaptive Feature Detection

```cpp
std::vector<Real> FindAllFeatures(const IRealFunction& func, Real a, Real b)
{
    RealFunctionAnalyzer analyzer(func);
    std::vector<Real> features;
    
    // Start with coarse sampling
    auto roots = analyzer.GetRoots(a, b, 1e-6);
    features.insert(features.end(), roots.begin(), roots.end());
    
    auto extrema = analyzer.GetLocalOptimums(a, b, 1e-6);
    features.insert(features.end(), extrema.begin(), extrema.end());
    
    auto inflections = analyzer.GetInflectionPoints(a, b, 1e-6);
    features.insert(features.end(), inflections.begin(), inflections.end());
    
    // Sort and remove duplicates
    std::sort(features.begin(), features.end());
    features.erase(std::unique(features.begin(), features.end(), 
                               [](Real a, Real b) { return std::abs(a - b) < 1e-6; }), 
                   features.end());
    
    return features;
}
```

---

## Summary

The Function Analyzer framework provides:

✅ **Critical Point Analysis** - Classify local min/max/saddle points with second derivative test  
✅ **Discontinuity Detection** - Find and classify JUMP, REMOVABLE, INFINITE, OSCILLATORY  
✅ **Limit Computation** - Numerical left/right limits with adaptive convergence  
✅ **Root Finding** - Comprehensive zero detection with bracketing + refinement  
✅ **Inflection Points** - Locate concavity changes via f''(x) sign changes  
✅ **Continuity Tests** - Point and interval verification with multi-stage analysis  
✅ **Function Comparison** - Absolute, relative, and integrated difference metrics  
✅ **Formatted Reports** - Print methods for analysis results  

**Related Documentation**:
- [Root_finding.md](Root_finding.md) - Underlying root finding algorithms
- [../core/Derivation.md](../core/Derivation.md) - Numerical differentiation used internally
- [../core/Integration.md](../core/Integration.md) - Integration for L² norms
- [Statistics.md](Statistics.md) - Statistical measures for periodic analysis
