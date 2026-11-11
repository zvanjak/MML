# Function Analyzer

**Numerical API**: `mml/algorithms/Analyzers/FunctionsAnalyzer.h`

**Optional reports**: `mml/tools/analyzers/FunctionAnalysisReports.h`
**Dependencies**: RootFinding, Derivation, Integration, Statistics

---

## Overview

The Function Analyzer framework provides comprehensive tools for analyzing real-valued functions: finding roots, extrema, inflection points, detecting discontinuities, and comparing functions. It combines numerical techniques from calculus, root finding, and statistics to provide deep insights into function behavior.

> **These are sampling heuristics, not proofs.** Every predicate and search evaluates the function on a
> fixed grid (`numPoints`) and probes it with finite-difference derivatives and absolute thresholds
> (`eps`). `IsMonotonicOnSamples` / `AppearsContinuous` returning `true` means *no counter-evidence was found at
> the sample points*; features narrower than the grid spacing are missed and steep-but-continuous
> functions can be flagged as discontinuous. Use the analyzer to get a first picture of an unfamiliar
> function - not as a correctness check inside an algorithm (CONTRACT_POLICY §11 / analysis v4 A13).

### Key Capabilities

- **Critical Point Analysis**: Locate and classify local minima, maxima, and saddle points
- **Discontinuity Detection**: Find and classify jump, removable, infinite, and oscillatory discontinuities
- **Continuity Exploration**: Probe for evidence of discontinuities at points and over intervals
- **Root Finding**: Locate zeros of functions in intervals
- **Inflection Points**: Find where concavity changes (second derivative sign changes)
- **Monotonicity Testing**: Check monotonic ordering of sampled values
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
| `RealFunctionAnalyzer` | Single function analysis | `GetRoots`, `GetLocalOptimums`, `FindDiscontinuityCandidates` |
| `RealFunctionComparer` | Compare two functions | `getAbsDiffMax`, `getIntegratedAbsDiff` |
| `ScalarFieldAnalyzer` | 3D scalar field analysis | `GradientAt`, `LaplacianAt`, `FindCriticalPoints`, `IsHarmonic` |
| `VectorFieldAnalyzer` | 3D vector field analysis | `DivergenceAt`, `CurlAt`, `IsSolenoidal`, `IsIrrotational`, `CirculationAroundCircle` |

### Discontinuity Types

| Type | Description | Example |
|------|-------------|---------|
| `JUMP` | Left/right limits differ | Step function, `floor(x)` |
| `REMOVABLE` | Limits match but f(x) differs/undefined | `sin(x)/x` at x=0 |
| `Infinite` | At least one limit is infinite | `1/x` at x=0 |
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
};
```

**Parameters:**
- `f` - Function implementing `IRealFunction` interface

Display names belong to report functions, not the numerical analyzer. Mathematical operations
return data without writing output; including the numerical header does not include reporting tools.

**Example 1: Basic Setup**

```cpp
#include <mml/algorithms/Analyzers/FunctionsAnalyzer.h>

// Polynomial: f(x) = x³ - 6x² + 9x + 1
auto polynomial = [](Real x) { return x*x*x - 6*x*x + 9*x + 1; };
RealFunctionFromStdFunc func(polynomial);

RealFunctionAnalyzer analyzer(func);
```

---

## Root Finding

### GetRoots - Find All Zeros in Interval

```cpp
std::vector<Real> GetRoots(Real x1, Real x2, Real eps = 1e-6);
```

**Algorithm**:
1. Sample 1000 subintervals, including both endpoints
2. Record samples whose absolute function value is below `eps`
3. Detect finite, nonzero opposite-sign neighbors without multiplying their values
4. Refine sign-change brackets with bisection to tolerance `eps`

This is a sampling search, not an all-roots guarantee: tangent roots and multiple roots
inside one grid cell can be missed. The tolerance is used both for sampled residuals and
for bisection, so rescaling the function also changes its residual-based stopping behavior.

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
RealFunctionAnalyzer analyzer(func);

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
bool IsDefinedAtPoint(Real x) const;
bool IsContinuousAtPoint(Real x, Real eps = 1e-6);
bool IsDerivativeDefinedAtPoint(Real x, Real eps = 1e-6);
bool IsInflectionPoint(Real x, Real eps = 1e-6);
bool IsLocalOptimum(Real x, Real eps = 1e-6);
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
std::cout << "  Defined: " << (analyzer.IsDefinedAtPoint(x) ? "yes" : "no") << std::endl;
std::cout << "  Continuous: " << (analyzer.IsContinuousAtPoint(x, 1e-6) ? "yes" : "no") << std::endl;
std::cout << "  Derivative exists: " << (analyzer.IsDerivativeDefinedAtPoint(x, 1e-6) ? "yes" : "no") << std::endl;

// Output:
// Defined: yes
// Continuous: yes
// Derivative exists: no (left deriv = -1, right deriv = 0)
```

### Interval-Level Tests

```cpp
bool AppearsContinuous(Real x1, Real x2, int numPoints, Real eps = 1e-6,
                  std::vector<Real>* discontinuities = nullptr);

bool IsMonotonicOnSamples(Real x1, Real x2, int numPoints);
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
bool continuous = analyzer.AppearsContinuous(0.0, 5.0, 500, 1e-6, &discontinuities);

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

Both directions preserve the sign of detected infinite limits. Finite values above
`Precision::InfiniteLimitThreshold` are not automatically infinite: the probe requires
repeated same-sign growth above that threshold, or repeated same-sign infinities.
A large finite constant therefore retains its finite estimate. These remain heuristic
estimates; exhausting the budget returns the last finite sample, not a convergence
certificate. NaN evaluations, unresolved nonfinite sequences, or the absence of a
representable one-sided probe return NaN. A probe never substitutes `f(x)` when its
step rounds back to `x`.

`ClassifyDiscontinuity` checks for infinite limits before missing estimates, including
when both sides are infinite. Missing/NaN estimates are `UNKNOWN`, not evidence of
oscillation; `OSCILLATORY` remains an available classification value but is not inferred
merely from invalid evaluations. Large finite estimates alone do not establish a pole.

The `SafeDifference` utility computes an absolute difference over the representable
range. It returns infinity for an infinite input or arithmetic overflow, and propagates
NaN inputs rather than inventing an infinite difference.

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

Use `DiscontinuityType::Infinite` for infinite limits. The numerical analyzer and optional
report header preserve caller-defined `INFINITE` and `DOMAIN` macros, including Windows
definitions. The deprecated `DiscontinuityType::INFINITE` alias retains the same value for
one minor release, but that legacy spelling is usable only when the `INFINITE` macro is
absent. No enumerator values or report labels have changed.

```cpp
DiscontinuityPoint ClassifyDiscontinuity(Real x, Real eps = 1e-6) const;
```

**DiscontinuityPoint Structure**:
```cpp
struct DiscontinuityPoint
{
    Real x;                       // Location
    DiscontinuityType type;       // JUMP, REMOVABLE, Infinite, OSCILLATORY, UNKNOWN
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
    case DiscontinuityType::Infinite:
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

std::cout << "Type: " << (disc.type == DiscontinuityType::Infinite ? "INFINITE" : "OTHER") << std::endl;
std::cout << "Left limit:  " << disc.leftLimit << std::endl;   // -∞
std::cout << "Right limit: " << disc.rightLimit << std::endl;  // +∞
```

### FindDiscontinuityCandidates - Comprehensive Search

```cpp
std::vector<DiscontinuityPoint> FindDiscontinuityCandidates(Real x1, Real x2,
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
auto discontinuities = analyzer.FindDiscontinuityCandidates(-Constants::PI, Constants::PI, 2000);

std::cout << "Found " << discontinuities.size() << " discontinuities:" << std::endl;
for (const auto& disc : discontinuities)
{
    std::cout << "\nAt x ≈ " << disc.x << ":" << std::endl;
    std::cout << "  Type: ";
    switch (disc.type)
    {
        case DiscontinuityType::Infinite:
            std::cout << "INFINITE (vertical asymptote)" << std::endl;
            break;
        default:
            std::cout << "OTHER" << std::endl;
    }
}
// Expected: discontinuities near -π/2, π/2
```

### AnalyzeContinuity and tools-side PrintContinuityAnalysis

`AnalyzeContinuity` returns candidate classifications together with the interval, subdivision
count, and probe tolerance. The report consumes this snapshot; it does not evaluate the function.

**Example 16: Full Continuity Report**

```cpp
#include <mml/tools/analyzers/FunctionAnalysisReports.h>

RealFunctionFromStdFunc f([](Real x) { return std::floor(x); });
const RealFunctionAnalyzer analyzer(f);
const auto continuity = analyzer.AnalyzeContinuity(-2.5, 2.5, 100);
AnalysisReports::PrintContinuityAnalysis(std::cout, continuity, "floor(x)");
```

An empty candidate list means **no candidates were detected**, not that continuity was proved.
The search uses `numSubintervals + 1` initial grid nodes and additional probes; the subdivision
count is not a total function-evaluation count. Invalid bounds, nonpositive subdivisions, and
nonfinite or nonpositive probe tolerances are rejected with `std::invalid_argument`.

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

## Optional Tools-Side Reports

All report functions live in `MML::AnalysisReports`, accept an explicit `std::ostream&`,
and consume already computed data. They never call the original function. They use classic-locale
six-significant-digit display formatting and restore the caller's locale, flags, precision, width,
and fill character, including when stream output throws. Standard ostream error-state and exception
mask behavior is preserved.

### PrintPointAnalysis

**Example 20: Point Report**

```cpp
auto func = [](Real x) { return std::abs(x); };
RealFunctionFromStdFunc f(func);
RealFunctionAnalyzer analyzer(f);

const auto point = analyzer.AnalyzePoint(0.0);
AnalysisReports::PrintPointAnalysis(std::cout, point, "|x|");
```

### PrintIntervalAnalysis

**Example 21: Interval Report**

```cpp
auto func = [](Real x) { return 1.0 / (1.0 + x*x); };  // Lorentzian
RealFunctionFromStdFunc f(func);
RealFunctionAnalyzer analyzer(f);

const auto interval = analyzer.AnalyzeInterval(-5.0, 5.0, 1000);
AnalysisReports::PrintIntervalAnalysis(std::cout, interval, "1/(1+x^2)");
```

The interval grid is unchanged: `x1 + i*(x2-x1)/numPoints`, for `0 <= i < numPoints`.
The right endpoint is excluded. `minValue` and `maxValue` are **finite sampled extrema**,
not global bounds. `AnalyzePoints` uses the same half-open grid. The structured interval
monotonicity flag is computed on that grid, not on a different endpoint set.
Any nonfinite sample makes the defined, continuity, and monotonicity flags false and
is listed among both undefined points and discontinuity candidates. If every sample is
nonfinite, both extrema are NaN and the report says extrema are unavailable rather than
printing numerical sentinel values.

In contrast, `IsMonotonicOnSamples`, `AppearsContinuous`, `MinInNPoints`,
`MaxInNPoints`, and `FindDiscontinuityCandidates` interpret their legacy `numPoints`
argument as a **subinterval count** and inspect its `numPoints + 1` grid nodes,
including both endpoints. Continuity searches may add refinement probes.
Monotonicity permits plateaus, determines direction from the first unequal pair, and
rejects nonfinite samples. Minimum/maximum searches ignore nonfinite samples, accept
every finite magnitude, and return NaN only when no finite sample exists.
`AppearsContinuous` replaces a supplied discontinuity vector on every successful
return, including clearing it when no candidates were found.

These sampling APIs require positive counts and finite ordered bounds (`x1 <= x2`);
tolerance-taking analysis/probe APIs require a positive finite `eps`.
Limit probes additionally require a positive iteration budget and a representable
initial step. Invalid parameters throw `std::invalid_argument` before sampling or
reserving output storage. A zero-width interval is permitted.

### Reusing Results and Detailed Tables

```cpp
RealFunction f(+[](Real x) { return x * x * x - 3 * x + 1; });
const RealFunctionAnalyzer analyzer(f);
const auto interval = analyzer.AnalyzeInterval(-2.0, 2.0, 100);
const auto details = analyzer.AnalyzePoints(-2.0, 2.0, 20);
AnalysisReports::PrintIntervalAnalysis(std::cout, interval);
AnalysisReports::PrintDetailedIntervalAnalysis(std::cout, details);
std::ostringstream savedReport;
AnalysisReports::PrintIntervalAnalysis(savedReport, interval);
```

The same snapshot can be rendered repeatedly without repeating any function evaluations.
For numerical point results, `isContinuous`, `isDerivativeDefined`, and `isInflectionPoint`
record probe outcomes, not existence proofs. Interval flags likewise describe sampled evidence;
their legacy field names do not strengthen these guarantees.

### Migration

- Pass a function only to `RealFunctionAnalyzer`; pass display names to report functions.
- Replace member `Print*` calls with `AnalysisReports::Print*(stream, result, optionalName)`.
- Use `IsMonotonicOnSamples`, `AppearsContinuous`, and `FindDiscontinuityCandidates`.
  The former `IsMonotonic`, `IsContinuous`, and `FindDiscontinuities` spellings remain deprecated
  numerical aliases. Printer member functions are removed, not forwarded back into tools.
- Numerical comparisons and searches remain in algorithms. Checked field sampling is documented below.

Runnable counterparts: [docs_demo_function_analyzer.cpp](../../src/docs_demos/algorithms/docs_demo_function_analyzer.cpp).
Build `MML_DocsApp` and run it with `--analyzer-reports`.

---

## Advanced Features

### ScalarFieldAnalyzer - 3D Field Analysis

```cpp
ScalarFunctionFromStdFunc<3> distSqField([](const VectorN<Real, 3>& p) {
    return p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
});
ScalarFieldAnalyzer scalarAnalyzer(distSqField);
Vec3Cart lower{-1, -1, -1}, upper{1, 1, 1};
const auto criticalPoints = scalarAnalyzer.FindCriticalPoints(lower, upper, 2);
const bool harmonic = scalarAnalyzer.IsHarmonic(lower, upper, 2);
std::cout << "Critical samples: " << criticalPoints.size()
          << ", harmonic: " << harmonic << '\n'; // 1, false
```

`GradientAt`, `IsCriticalPoint`, and `LaplacianAt` evaluate individual positions.
`FindCriticalPoints` returns grid points whose gradient magnitude is below tolerance;
it does not refine roots or prove that other critical points are absent.
`IsHarmonic` checks the absolute Laplacian against tolerance at each visited grid point.

### VectorFieldAnalyzer - 3D Vector Analysis

```cpp
VectorFunctionFromStdFunc<3> rotationField([](const VectorN<Real, 3>& p) {
    return VectorN<Real, 3>{p[1], -p[0], 0};
});
VectorFieldAnalyzer vectorAnalyzer(rotationField);
Vec3Cart lower{-1, -1, -1}, upper{1, 1, 1};
const bool solenoidal = vectorAnalyzer.IsSolenoidal(lower, upper, 2);
const bool curlFreeOnSamples = vectorAnalyzer.IsConservative(lower, upper, 2);
const Real maxCurl = vectorAnalyzer.MaxCurlMagnitude(lower, upper, 2);
std::cout << "Solenoidal: " << solenoidal
          << ", curl-free on samples: " << curlFreeOnSamples
          << ", maximum curl magnitude: " << maxCurl << '\n'; // true, false, about 2
```

`DivergenceAt`, `CurlAt`, `IsSolenoidalAt`, and `IsIrrotationalAt` evaluate a point.
`MaxDivergence` and `MaxCurlMagnitude` return sampled maxima, not certified bounds.
`CirculationAroundCircle` approximates the oriented circular line integral using
midpoint field values on polygonal chords.

### Sampling, failure, and lifetime contracts

- Both analyzers **borrow** a const field reference. The field must outlive the analyzer;
  constructing an analyzer directly from a temporary field is rejected at compile time.
- `samplesPerDim` is a positive number of subintervals. Each axis includes both endpoints:
  a completed traversal visits `(samplesPerDim + 1)^3` nodes. Equal endpoints are allowed
  and retain repeated samples on collapsed axes; reversed bounds are invalid.
  Coordinates are interpolated without forming a potentially overflowing endpoint difference.
- Coordinates and bounds must be finite, and every explicit tolerance must be positive
  and finite. Invalid numeric inputs throw `NumericInputError`; nonpositive counts or
  reversed bounds throw `ArgumentError`. An unrepresentable grid node count throws
  `ResourceLimitError`, rather than overflowing loop or work accounting.
- Differential evaluations use the checked Cartesian
  [Core field operations](../core/Field_operations.md), preserving default stencils and
  derivative-provider dispatch without requesting error estimates. Nonfinite consumed
  samples, provider derivatives, resulting derivatives, or magnitude/arithmetic reductions
  throw `NumericalMethodError`; they never become successful flags or zero maxima.
  Exceptions thrown by field callbacks propagate unchanged.
- Region predicates short-circuit on the first failing sample. A successful result
  describes only the sampled evidence; a negative result need not evaluate later samples.
- Circle centers/normals must be finite, normals nonzero, radii finite and nonnegative,
  and `numPoints` positive. Invalid numeric inputs throw `NumericInputError`, while
  negative radii, zero normals, or nonpositive counts throw `ArgumentError`.
  A zero radius is supported. Finite samples and representable path/integral arithmetic
  are required even for a degenerate circle.
- `IsConservative` is retained as a compatibility spelling for **sampled irrotationality**
  (`IsIrrotational`). It does not establish global path independence or even zero curl
  between samples. The mathematical implication from curl-free to conservative requires
  appropriate smoothness and domain topology (for example a simply connected domain).

Runnable counterparts: [docs_demo_field_analyzers.cpp](../../src/docs_demos/algorithms/docs_demo_field_analyzers.cpp).
Build `MML_DocsApp` and run it with `--field-analyzers`.

---

## Integration with Other Modules

### With Root Finding

```cpp
#include <mml/algorithms/RootFinding.h>
#include <mml/algorithms/Analyzers/FunctionsAnalyzer.h>

// Function analyzer uses RootFinding internally
RealFunctionAnalyzer analyzer(func);
std::vector<Real> roots = analyzer.GetRoots(a, b);  // Uses RootFinding::FindRootBrackets + Bisection
```

### With Derivation

```cpp
#include <mml/core/Derivation.h>

// All derivative-based methods use Derivation module
Real firstDeriv = Derivation::NDer4(func, x);
Real secondDeriv = Derivation::NSecDer4(func, x);

// Used internally for critical points and inflection detection
```

### With Integration

```cpp
#include <mml/core/Integration.h>

// Function comparer uses integration for L² norms
RealFunctionComparer comparer(f1, f2);
Real sqrDiff = comparer.getIntegratedSqrDiff(a, b, IntegrationMethod::SIMPSON);
```

### With Statistics

```cpp
#include <mml/algorithms/Statistics.h>

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
1. `FindDiscontinuityCandidates` - multi-stage refinement
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
    RealFunctionAnalyzer analyzer(func);
    
    // Find all interesting features
    auto roots = analyzer.GetRoots(a, b);
    auto extrema = analyzer.GetLocalOptimumsClassified(a, b);
    auto inflections = analyzer.GetInflectionPoints(a, b);
    auto discontinuities = analyzer.FindDiscontinuityCandidates(a, b);
    
    // Print comprehensive report
    std::cout << "=== FUNCTION ANALYSIS ===" << std::endl;
    std::cout << "Roots: " << roots.size() << std::endl;
    std::cout << "Local extrema: " << extrema.size() << std::endl;
    std::cout << "Inflection points: " << inflections.size() << std::endl;
    std::cout << "Discontinuities: " << discontinuities.size() << std::endl;
    
    AnalysisReports::PrintContinuityAnalysis(std::cout, analyzer.AnalyzeContinuity(a, b));
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
✅ **Discontinuity Detection** - Find and classify JUMP, REMOVABLE, Infinite, OSCILLATORY
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
