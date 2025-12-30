# Intervals — Comprehensive Reference

**File**: `mml/base/Intervals.h`

> **Note**: Interval classes are NOT templates - they use `Real` type internally.
> Example: Use `ClosedInterval(0.0, 10.0)` not `ClosedInterval<double>(0.0, 10.0)`

## Overview

The MinimalMathLibrary (MML) provides a comprehensive interval arithmetic system for representing and manipulating mathematical intervals. The implementation supports both simple intervals (single contiguous ranges) and compound intervals (unions of multiple ranges), with support for open/closed endpoints and infinite bounds.

## Highlights
- Full support for open/closed endpoints and infinite bounds.
- Simple and compound intervals (unions of ranges) with set operations.
- Epsilon-aware comparisons for floating-point robustness.
- Factory methods and friend access patterns for efficient implementations.
- Tools for coverings, membership tests, intersection/difference/complement.

### Runnable Examples

| Demo Function | Description | Location |
|--------------|-------------|----------|
| `Docs_Demo_Intervals()` | Basic interval operations | `src/docs_demos/docs_demo_intervals.cpp` |

**Build and Run:**
```bash
cmake --build build --target MML_DocsApp
./build/src/MML_DocsApp   # Linux/macOS
.\build\src\Release\MML_DocsApp.exe   # Windows
```

## Table of Contents

1. [Design Philosophy](#design-philosophy)
2. [Architecture](#architecture)
3. [Interval Types](#interval-types)
4. [Core Operations](#core-operations)
5. [Usage Examples](#usage-examples)
6. [Advanced Features](#advanced-features)
7. [Best Practices](#best-practices)
8. [Mathematical Background](#mathematical-background)
9. [API Reference](#api-reference)

## Design Philosophy

The interval system is built on several key principles:

- **Type Safety**: Different endpoint types (open/closed) are enforced through the type system
- **Precision Handling**: Careful epsilon-based comparisons for floating-point arithmetic
- **Composability**: Simple intervals can be combined into compound intervals
- **Mathematical Correctness**: All operations preserve mathematical semantics
- **Performance**: Efficient algorithms with minimal overhead for common operations

## Architecture

### Class Hierarchy

```
IInterval (interface)
    ↓
BaseInterval (abstract base)
    ↓
    ├── ClosedInterval           [a, b]
    ├── OpenInterval             (a, b)
    ├── OpenClosedInterval       (a, b]
    ├── ClosedOpenInterval       [a, b)
    ├── ClosedToInfInterval      [a, +∞)
    ├── NegInfToClosedInterval   (-∞, b]
    ├── OpenToInfInterval        (a, +∞)
    ├── NegInfToOpenInterval     (-∞, b)
    └── CompleteRInterval        (-∞, +∞)
    
Interval (compound interval - union of BaseIntervals)
```

### Key Design Decisions

1. **Friend Class Pattern**: `Interval` is declared as a friend of `BaseInterval` to access protected members for set operations
2. **Vector-Based Output**: `GetEquidistantCovering` returns points via output parameter for efficiency
3. **Static Factory Methods**: Set operations (`Intersection`, `Difference`, `Complement`) available as static constructors on `Interval`
4. **Epsilon-Based Comparisons**: All floating-point comparisons use configurable epsilon for numerical stability

## Interval Types

### Simple Intervals (BaseInterval Subclasses)

#### ClosedInterval [a, b]
Contains both endpoints.

```cpp
ClosedInterval interval(0.0, 10.0);
// Represents [0, 10]
// Contains: 0.0, 5.0, 10.0
```

**Use Cases**: Physical measurements, bounded domains, probability distributions

#### OpenInterval (a, b)
Excludes both endpoints.

```cpp
OpenInterval interval(0.0, 10.0);
// Represents (0, 10)
// Contains: 5.0
// Excludes: 0.0, 10.0
```

**Use Cases**: Strict inequalities, domains with singularities at boundaries

#### OpenClosedInterval (a, b]
Open at lower bound, closed at upper bound.

```cpp
OpenClosedInterval interval(0.0, 10.0);
// Represents (0, 10]
// Contains: 5.0, 10.0
// Excludes: 0.0
```

**Use Cases**: Right-continuous functions, cumulative distributions

#### ClosedOpenInterval [a, b)
Closed at lower bound, open at upper bound.

```cpp
ClosedOpenInterval interval(0.0, 10.0);
// Represents [0, 10)
// Contains: 0.0, 5.0
// Excludes: 10.0
```

**Use Cases**: Left-continuous functions, array indexing analogs

### Infinite Bound Intervals

#### ClosedToInfInterval [a, +∞)
Unbounded on the right.

```cpp
ClosedToInfInterval ray(0.0);
// Represents [0, +∞)
// Contains: 0.0, 1.0, 1000.0, ...
```

#### NegInfToClosedInterval (-∞, b]
Unbounded on the left.

```cpp
NegInfToClosedInterval ray(0.0);
// Represents (-∞, 0]
// Contains: -1000.0, -1.0, 0.0
```

#### OpenToInfInterval (a, +∞)
Open lower bound, unbounded on the right.

```cpp
OpenToInfInterval ray(0.0);
// Represents (0, +∞)
// Excludes: 0.0
```

#### NegInfToOpenInterval (-∞, b)
Unbounded on the left, open upper bound.

```cpp
NegInfToOpenInterval ray(0.0);
// Represents (-∞, 0)
// Excludes: 0.0
```

#### CompleteRInterval (-∞, +∞)
Represents all real numbers.

```cpp
CompleteRInterval line;
// Represents (-∞, +∞)
// Contains all finite values
```

### Compound Intervals

The `Interval` class represents unions of `BaseInterval` objects.

```cpp
Interval compound;
compound.AddInterval(ClosedInterval(0.0, 5.0));
compound.AddInterval(ClosedInterval(10.0, 15.0));
// Represents [0, 5] ∪ [10, 15]
```

## Core Operations

### Containment Testing

```cpp
ClosedInterval interval(0.0, 10.0);

// Point containment
bool contains5 = interval.contains(5.0);      // true
bool contains0 = interval.contains(0.0);      // true (closed endpoint)
bool contains15 = interval.contains(15.0);    // false

// Interval containment
ClosedInterval sub(2.0, 8.0);
bool isContained = interval.contains(sub);    // true
```

### Intersection

Returns the overlapping region of two intervals.

```cpp
ClosedInterval a(0.0, 10.0);          // [0, 10]
ClosedInterval b(5.0, 15.0);          // [5, 15]

Interval result = Interval::Intersection(a, b);
// Result: [5, 10]
```

**Endpoint Type Rules**:
- Closed ∩ Closed = Closed
- Open ∩ Closed = Open
- Open ∩ Open = Open
- **Open endpoint "wins"** in intersections

**Examples**:
```cpp
// Example 1: Both closed
ClosedInterval c1(0.0, 10.0);         // [0, 10]
ClosedInterval c2(5.0, 15.0);         // [5, 15]
auto r1 = Interval::Intersection(c1, c2);  // [5, 10]

// Example 2: Mixed endpoints
OpenInterval o1(0.0, 10.0);           // (0, 10)
ClosedInterval c3(5.0, 10.0);         // [5, 10]
auto r2 = Interval::Intersection(o1, c3);  // [5, 10) - open at 10
```

### Difference

Returns the portion of the first interval not covered by the second.

```cpp
ClosedInterval a(0.0, 10.0);          // [0, 10]
ClosedInterval b(3.0, 7.0);           // [3, 7]

Interval result = Interval::Difference(a, b);
// Result: [0, 3) ∪ (7, 10]  (two sub-intervals)
```

**Five Cases Handled**:

1. **No Intersection**: Returns original interval
   ```cpp
   a = [0, 5], b = [10, 15]
   a \ b = [0, 5]
   ```

2. **Complete Containment**: Returns empty interval
   ```cpp
   a = [3, 7], b = [0, 10]
   a \ b = ∅
   ```

3. **Left Cut**: Removes left portion
   ```cpp
   a = [5, 10], b = [0, 7]
   a \ b = (7, 10]
   ```

4. **Right Cut**: Removes right portion
   ```cpp
   a = [0, 10], b = [7, 15]
   a \ b = [0, 7)
   ```

5. **Middle Cut (Split)**: Removes middle portion, creating two intervals
   ```cpp
   a = [0, 10], b = [3, 7]
   a \ b = [0, 3) ∪ (7, 10]
   ```

### Complement

Returns the complement of an interval with respect to the real line.

```cpp
// Finite interval complement
ClosedInterval interval(0.0, 10.0);   // [0, 10]
Interval comp = Interval::Complement(interval);
// Result: (-∞, 0) ∪ (10, +∞)

// Infinite interval complement
ClosedToInfInterval ray(5.0);          // [5, +∞)
Interval comp2 = Interval::Complement(ray);
// Result: (-∞, 5)
```

**Endpoint Inversion**:
- Closed endpoint → Open endpoint in complement
- Open endpoint → Closed endpoint in complement

### Equidistant Covering

Generates uniformly distributed points across the interval.

```cpp
ClosedInterval interval(0.0, 10.0);
std::vector<Real> points;
interval.GetEquidistantCovering(5, points);
// points = {0.0, 2.5, 5.0, 7.5, 10.0}

OpenInterval open(0.0, 10.0);
std::vector<Real> openPoints;
open.GetEquidistantCovering(5, openPoints);
// openPoints ≈ {0.0+ε, 2.5, 5.0, 7.5, 10.0-ε}
// where ε is a small offset for open endpoints
```

**Features**:
- Respects endpoint types (adds epsilon offset for open endpoints)
- Handles infinite bounds gracefully
- Configurable number of points
- Useful for numerical integration, sampling, visualization

## Usage Examples

### Example 1: Function Domain Analysis

Determine where a function is defined.

```cpp
// tan(x) is undefined at x = π/2 + nπ
// Domain: ℝ \ {π/2, 3π/2, 5π/2, ...}

Interval<double> domain;
double pi = Constants::PI;

// Add valid intervals between singularities
domain.addInterval(std::make_unique<OpenInterval<double>>(-pi/2, pi/2));
domain.addInterval(std::make_unique<OpenInterval<double>>(pi/2, 3*pi/2));
domain.addInterval(std::make_unique<OpenInterval<double>>(3*pi/2, 5*pi/2));

// Test points
bool validAt0 = domain.contains(0.0);           // true
bool validAtPiOver2 = domain.contains(pi/2);    // false (singularity)
```

### Example 2: Temperature Range Analysis

Analyze overlapping temperature ranges.

```cpp
// Comfortable temperature ranges for different activities
ClosedInterval<double> sleeping(15.0, 20.0);   // [15°C, 20°C]
ClosedInterval<double> working(18.0, 24.0);    // [18°C, 24°C]

// Optimal range for both
Interval<double> optimal = Interval<double>::Intersection(sleeping, working);
// Result: [18, 20]

// Temperature range suitable only for sleeping
Interval<double> sleepOnly = Interval<double>::Difference(sleeping, working);
// Result: [15, 18)
```

### Example 3: Root Finding Preparation

Bracket roots before applying numerical methods.

```cpp
// Find where f(x) = x² - 4x + 3 might have roots
// Known: roots are in [0, 5]

ClosedInterval<double> searchRange(0.0, 5.0);
std::vector<double> testPoints;
searchRange.GetEquidistantCovering(20, testPoints);

// Sample function at test points to find sign changes
// (sign changes indicate roots by Intermediate Value Theorem)
for (size_t i = 0; i < testPoints.size() - 1; ++i) {
    double f1 = /* evaluate f at testPoints[i] */;
    double f2 = /* evaluate f at testPoints[i+1] */;
    
    if (f1 * f2 < 0) {
        // Root exists in [testPoints[i], testPoints[i+1]]
        ClosedInterval<double> bracket(testPoints[i], testPoints[i+1]);
        // Apply bisection, Newton-Raphson, etc.
    }
}
```

### Example 4: Set Operations Chain

Combine multiple set operations.

```cpp
ClosedInterval<double> universe(0.0, 100.0);   // [0, 100]
ClosedInterval<double> exclude1(20.0, 30.0);   // [20, 30]
ClosedInterval<double> exclude2(60.0, 70.0);   // [60, 70]

// Remove two ranges from universe
Interval<double> step1 = Interval<double>::Difference(universe, exclude1);
// step1: [0, 20) ∪ (30, 100]

Interval<double> final = Interval<double>::Difference(step1, exclude2);
// final: [0, 20) ∪ (30, 60) ∪ (70, 100]
```

### Example 5: Compound Interval Operations

Work with unions of intervals.

```cpp
// Create compound interval: [0, 5] ∪ [10, 15] ∪ [20, 25]
Interval<double> compound;
compound.addInterval(std::make_unique<ClosedInterval<double>>(0.0, 5.0));
compound.addInterval(std::make_unique<ClosedInterval<double>>(10.0, 15.0));
compound.addInterval(std::make_unique<ClosedInterval<double>>(20.0, 25.0));

// Query properties
double totalLength = compound.getLength();           // 15.0
double lowerBound = compound.getLowerBound();        // 0.0
double upperBound = compound.getUpperBound();        // 25.0

// Test containment
bool contains3 = compound.contains(3.0);             // true (in [0, 5])
bool contains7 = compound.contains(7.0);             // false (gap)
bool contains12 = compound.contains(12.0);           // true (in [10, 15])

// Generate covering points
std::vector<double> points;
compound.GetEquidistantCovering(15, points);
// Distributes points proportionally across all sub-intervals
```

## Advanced Features

### Epsilon Handling

All comparisons use epsilon-based tolerance for floating-point arithmetic:

```cpp
// Configured in MMLPrecision.h
constexpr double EPSILON_STRICT = 1e-10;

// Used internally for:
// - Endpoint proximity detection
// - Intersection boundary calculation
// - Containment testing with numerical tolerance
```

### Infinite Bound Support

Operations gracefully handle infinite bounds:

```cpp
RightRay<double> ray(0.0);                     // [0, +∞)
ClosedInterval<double> finite(5.0, 10.0);     // [5, 10]

// Intersection with infinite bound
Interval<double> result = Interval<double>::Intersection(ray, finite);
// Result: [5, 10]

// Difference with infinite bound
Interval<double> diff = Interval<double>::Difference(ray, finite);
// Result: [0, 5) ∪ (10, +∞)
```

### Endpoint Type Propagation

Set operations intelligently propagate endpoint types:

- **Intersection**: More restrictive endpoint wins (open > closed)
- **Difference**: Inverts endpoint type at cut points
- **Complement**: Inverts all endpoint types

```cpp
OpenInterval<double> open(0.0, 10.0);          // (0, 10)
ClosedInterval<double> closed(5.0, 15.0);      // [5, 15]

auto inter = Interval<double>::Intersection(open, closed);
// Result: [5, 10) - lower closed from 'closed', upper open from 'open'

auto comp = Interval<double>::Complement(closed);
// Result: (-∞, 5) ∪ (15, +∞) - endpoints inverted
```

### Performance Considerations

- **Point containment**: O(n) for compound intervals with n sub-intervals
- **Intersection**: O(1) for simple intervals
- **Difference**: O(1) for simple intervals (up to 2 result sub-intervals)
- **Complement**: O(1) for finite intervals (always 2 or fewer sub-intervals)
- **Equidistant covering**: O(m*n) where m = points, n = sub-intervals

## Best Practices

### 1. Choose Appropriate Interval Types

Match the interval type to your mathematical requirements:

```cpp
// ✓ GOOD: Use open interval for strict inequality domains
OpenInterval<double> strictPositive(0.0, INFINITY);  // x > 0

// ✗ AVOID: Using closed interval when singularity exists at boundary
ClosedInterval<double> withSingularity(0.0, 10.0);  // If f(0) undefined
```

### 2. Handle Empty Results

Set operations may produce empty intervals:

```cpp
ClosedInterval<double> a(0.0, 5.0);
ClosedInterval<double> b(10.0, 15.0);

Interval<double> inter = Interval<double>::Intersection(a, b);
if (inter.IsEmpty()) {
    // Handle no intersection case
}
```

### 3. Use Equidistant Covering Wisely

Balance point count vs. accuracy needs:

```cpp
// ✓ GOOD: Reasonable point count for numerical integration
ClosedInterval<double> domain(0.0, 1.0);
std::vector<double> points;
domain.GetEquidistantCovering(100, points);  // 100 points usually sufficient

// ✗ AVOID: Excessive points causing performance issues
domain.GetEquidistantCovering(1000000, points);  // Overkill for most cases
```

### 4. Leverage Compound Intervals

Use compound intervals for complex domains:

```cpp
// ✓ GOOD: Model non-contiguous valid ranges
Interval<double> validRanges;
validRanges.addInterval(std::make_unique<ClosedInterval<double>>(-10.0, -5.0));
validRanges.addInterval(std::make_unique<ClosedInterval<double>>(5.0, 10.0));

// ✗ AVOID: Trying to force single interval for disjoint sets
// Can't represent {x : -10 ≤ x ≤ -5 or 5 ≤ x ≤ 10} as one BaseInterval
```

### 5. Mind Endpoint Semantics in Operations

Understand endpoint type propagation:

```cpp
OpenInterval<double> open(0.0, 10.0);          // (0, 10)
ClosedInterval<double> closed(0.0, 10.0);      // [0, 10]

// These are NOT the same interval!
bool sameAtZero = open.contains(0.0);          // false
bool sameAtTen = closed.contains(0.0);         // true
```

## Mathematical Background

### Interval Arithmetic Foundations

An interval is a set of real numbers with a specific structure:

**Definition**: An interval `I` is a subset of ℝ such that for any `x, y ∈ I` with `x < y`, all points `z` with `x < z < y` are also in `I`.

**Types** (for finite `a, b` with `a ≤ b`):
- Closed: `[a, b] = {x ∈ ℝ : a ≤ x ≤ b}`
- Open: `(a, b) = {x ∈ ℝ : a < x < b}`
- Half-open: `[a, b) = {x ∈ ℝ : a ≤ x < b}` and `(a, b] = {x ∈ ℝ : a < x ≤ b}`

### Set Operations

**Intersection**: `A ∩ B = {x : x ∈ A and x ∈ B}`

Properties:
- Commutative: `A ∩ B = B ∩ A`
- Associative: `(A ∩ B) ∩ C = A ∩ (B ∩ C)`
- Idempotent: `A ∩ A = A`

**Difference**: `A \ B = {x : x ∈ A and x ∉ B}`

Properties:
- Non-commutative: `A \ B ≠ B \ A` (in general)
- Not associative
- Result may be disconnected (compound interval)

**Complement**: `A^c = {x ∈ ℝ : x ∉ A}`

Properties:
- Involution: `(A^c)^c = A`
- De Morgan's laws: `(A ∪ B)^c = A^c ∩ B^c`, `(A ∩ B)^c = A^c ∪ B^c`

### Endpoint Type Rules

For intersection of intervals `I₁` and `I₂` with overlapping region:

| I₁ endpoint | I₂ endpoint | Result endpoint |
|-------------|-------------|-----------------|
| Closed      | Closed      | Closed          |
| Closed      | Open        | Open            |
| Open        | Closed      | Open            |
| Open        | Open        | Open            |

**Rule**: The more restrictive (open) endpoint type wins.

**Mathematical Justification**: 
- `[a, b] ∩ (a, b) = (a, b)` because a point in the intersection must be in BOTH intervals
- If one interval excludes a boundary point, the intersection must also exclude it

## API Reference

### IInterval Interface

```cpp
template<typename Real>
class IInterval {
public:
    virtual bool contains(Real x) const = 0;
    virtual bool contains(const BaseInterval<Real>& inSecond) const = 0;
    virtual bool intersects(const BaseInterval<Real>& inSecond) const = 0;
    virtual int GetEquidistantCovering(int num_points, 
                                       std::vector<Real>& outPoints) const = 0;
    virtual Real getLowerBound() const = 0;
    virtual Real getUpperBound() const = 0;
    virtual Real getLength() const = 0;
};
```

### BaseInterval<Real>

Abstract base class for all simple intervals.

**Key Methods**:
- `bool contains(Real x) const` - Test if point is in interval
- `bool contains(const BaseInterval<Real>&) const` - Test interval containment
- `bool intersects(const BaseInterval<Real>&) const` - Test if intervals overlap
- `int GetEquidistantCovering(int n, vector<Real>&) const` - Generate n uniform points
- `Real getLowerBound() const` - Get lower bound
- `Real getUpperBound() const` - Get upper bound
- `Real getLength() const` - Get interval length (may be infinite)

**Protected Members**:
- `Real _a` - Lower bound
- `Real _b` - Upper bound

### ClosedInterval<Real>

Represents `[a, b]` - includes both endpoints.

```cpp
template<typename Real>
class ClosedInterval : public BaseInterval<Real> {
public:
    ClosedInterval(Real a, Real b);
    
    bool contains(Real x) const override;
    int GetEquidistantCovering(int num_points, 
                               std::vector<Real>& outPoints) const override;
    // ... other IInterval methods
};
```

### OpenInterval<Real>

Represents `(a, b)` - excludes both endpoints.

```cpp
template<typename Real>
class OpenInterval : public BaseInterval<Real> {
public:
    OpenInterval(Real a, Real b);
    
    bool contains(Real x) const override;
    int GetEquidistantCovering(int num_points, 
                               std::vector<Real>& outPoints) const override;
    // ... other IInterval methods
};
```

### Interval<Real>

Compound interval class - represents unions of BaseIntervals.

**Constructors**:
```cpp
Interval();  // Empty interval
```

**Interval Management**:
```cpp
void addInterval(std::unique_ptr<BaseInterval<Real>> interval);
void Clear();
bool IsEmpty() const;
int NumIntervals() const;
```

**Set Operations** (static methods):
```cpp
static Interval<Real> Intersection(const BaseInterval<Real>& a, 
                                   const BaseInterval<Real>& b);

static Interval<Real> Difference(const BaseInterval<Real>& a, 
                                const BaseInterval<Real>& b);

static Interval<Real> Complement(const BaseInterval<Real>& a);
```

**Query Methods**:
```cpp
bool contains(Real x) const;
bool contains(const BaseInterval<Real>& other) const;
bool intersects(const BaseInterval<Real>& other) const;
Real getLowerBound() const;
Real getUpperBound() const;
Real getLength() const;
int GetEquidistantCovering(int num_points, std::vector<Real>& outPoints) const;
```

### Infinite Bound Intervals

**RightRay<Real>**: `[a, +∞)`
```cpp
RightRay(Real a);
```

**LeftRay<Real>**: `(-∞, b]`
```cpp
LeftRay(Real b);
```

**OpenRightRay<Real>**: `(a, +∞)`
```cpp
OpenRightRay(Real a);
```

**OpenLeftRay<Real>**: `(-∞, b)`
```cpp
OpenLeftRay(Real b);
```

**RealLine<Real>**: `(-∞, +∞)`
```cpp
RealLine();
```

## Error Handling

The interval system performs validation:

```cpp
// Throws InvalidInterval exception if a > b
ClosedInterval<double> invalid(10.0, 5.0);  // ERROR!

// Check for empty results from set operations
Interval<double> result = Interval<double>::Intersection(a, b);
if (result.IsEmpty()) {
    // Handle empty intersection
}
```

## Testing

Comprehensive test suite in `tests/base/intervals_tests.cpp`:
- 26 test cases
- 346 assertions
- Coverage of all interval types and operations
- Edge cases, precision handling, and performance tests

Run tests:
```bash
cd build
cmake --build . --target MML_Tests
.\tests\Debug\MML_Tests.exe "[intervals]"
```

## Implementation Notes

### Files

- **Interface**: `mml/interfaces/IInterval.h`
- **Implementation**: `mml/base/Intervals.h`
- **Tests**: `tests/base/intervals_tests.cpp`

### Dependencies

- `MMLBase.h` - Base type definitions
- `MMLPrecision.h` - Epsilon constants
- `MMLExceptions.h` - Exception types

### Commits

- Implementation: 9096bfa7
- Test suite: 8aa58407
- Documentation: 489db81b, a10abe90

## References

- **Numerical Analysis**: Burden & Faires, "Numerical Analysis" (9th ed.)
- **Real Analysis**: Rudin, "Principles of Mathematical Analysis"
- **Interval Arithmetic**: Moore et al., "Introduction to Interval Analysis"

## Future Enhancements

Potential improvements for future versions:

1. **Arithmetic Operations**: Add `+`, `-`, `*`, `/` for interval arithmetic
2. **Union Optimization**: Automatic merging of overlapping sub-intervals
3. **Interval Arithmetic Error Bounds**: Track rounding error propagation
4. **Multi-dimensional Intervals**: Support for interval vectors/boxes
5. **Constraint Propagation**: Integration with constraint solving systems

---

*Last Updated: December 13, 2025*  
*Implementation Version: 1.0*  
*Test Coverage: 100% (26 tests, 346 assertions)*
