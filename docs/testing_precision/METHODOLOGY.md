# Precision Testing Methodology

## Overview

This document describes the methodology, framework design, and best practices used in MML's precision testing suite. It provides guidance for extending tests and interpreting results.

## Philosophy

### Why Precision Testing?

Numerical algorithms are deceptively fragile. An implementation can:
- Compile and run without errors
- Produce plausible-looking results
- Pass basic sanity checks
- **Yet still be wrong by orders of magnitude**

Precision testing ensures algorithms achieve their theoretical accuracy limits.

### Goals

1. **Verify convergence order** - Does NDer4 really give O(h⁴)?
2. **Quantify error bounds** - What accuracy can users expect?
3. **Compare algorithms** - Which method is best for which problem?
4. **Detect regressions** - Catch precision degradation during refactoring
5. **Document behavior** - Provide reference data for users

## The PrecisionTestFramework

### Architecture

```
PrecisionTestSuite
├── addResult(algorithm, function, exact, computed, time)
├── printSummary()
├── printErrorOrderMatrix()
├── exportCSV(filename)
└── getResults() -> vector<PrecisionTestResult>

PrecisionTestResult
├── algorithmName: string
├── functionName: string
├── exactValue: double
├── computedValue: double
├── absoluteError: double
├── relativeError: double
├── timingMs: double
└── errorOrder: int  // floor(log10(|error|))
```

### Usage Pattern

```cpp
#include "core/PrecisionTestFramework.h"

void Test_Precision_MyAlgorithm(MML::Tests::TestRegister& tests) {
    auto& suite = MML::Tests::GetPrecisionSuite();
    suite.clearResults();  // Start fresh
    
    // Test case 1
    double exact = 1.0;
    auto start = high_resolution_clock::now();
    double computed = myAlgorithm(params);
    auto end = high_resolution_clock::now();
    double time_ms = duration<double, milli>(end - start).count();
    
    suite.addResult("MyAlgorithm", "TestCase1", exact, computed, time_ms);
    
    // More test cases...
    
    // Output results
    suite.printSummary();
    suite.printErrorOrderMatrix();
    suite.exportCSV("results/precision_myalgorithm.csv");
}
```

## Test Design Principles

### 1. Use Known Exact Solutions

**Good:**
```cpp
// Derivative of sin(x) at x=1 is cos(1)
double exact = std::cos(1.0);
double computed = derivative(sin, 1.0, h);
```

**Bad:**
```cpp
// No known exact value
double computed = derivative(complex_function, 1.0, h);
// How do we validate?
```

### 2. Test Across Parameter Ranges

```cpp
// Test multiple step sizes
for (double h : {1e-2, 1e-4, 1e-6, 1e-8, 1e-10}) {
    double computed = derivative(f, x, h);
    suite.addResult("CentralDiff", "h=" + to_string(h), exact, computed, time);
}
```

### 3. Include Edge Cases

- Very small inputs (underflow risk)
- Very large inputs (overflow risk)
- Near singularities
- Boundary conditions
- Ill-conditioned problems

### 4. Test Convergence Order

For method with order p, verify:
```
error(h/2) / error(h) ≈ 2^(-p)
```

```cpp
double e1 = std::abs(compute(h) - exact);
double e2 = std::abs(compute(h/2) - exact);
double order = std::log2(e1 / e2);
// Should be approximately p
```

### 5. Compare Multiple Algorithms

```cpp
suite.addResult("Bisection", "cos(x)-x", exact, bisection(f, a, b), t1);
suite.addResult("Newton", "cos(x)-x", exact, newton(f, df, x0), t2);
suite.addResult("Brent", "cos(x)-x", exact, brent(f, a, b), t3);
// Now we can compare accuracy AND speed
```

## Error Metrics

### Absolute Error
```
ε_abs = |exact - computed|
```
Use when exact value is near 1 or when units matter.

### Relative Error
```
ε_rel = |exact - computed| / |exact|
```
Use when scale varies. **Undefined when exact = 0!**

### Error Order
```
order = floor(log10(|error|))
```
Quick visual indicator: -15 means ~10⁻¹⁵ accuracy.

### Residual (for linear systems)
```
residual = ||Ax - b||
```
Measures how well solution satisfies equations.

### Reconstruction Error (for decompositions)
```
||A - QR||, ||A - LDL^T||, ||A - UΣV^T||
```
Verifies decomposition quality.

## Standard Test Functions

### Derivatives
| Function | f(x) | f'(x) | Notes |
|----------|------|-------|-------|
| Polynomial | x³ - 2x² + x | 3x² - 4x + 1 | Exact arithmetic |
| Trigonometric | sin(x) | cos(x) | Well-behaved |
| Exponential | e^x | e^x | Self-derivative |
| Logarithmic | ln(x) | 1/x | Singularity at 0 |

### Integrals
| Function | ∫f dx | Notes |
|----------|-------|-------|
| Polynomial | ∫x² dx = x³/3 | Exact for high-order quadrature |
| Gaussian | ∫e^(-x²) dx = √π·erf | Standard test |
| Oscillatory | ∫sin(kx) dx | Tests sampling |
| Singular | ∫1/√x dx | Endpoint singularity |

### ODEs
| System | Exact Solution | Notes |
|--------|---------------|-------|
| y' = -y | y = e^(-t) | Exponential decay |
| y'' = -y | y = cos(t) | Harmonic oscillator |
| Lorenz | N/A | Chaotic, test structure |

### Linear Systems
| Matrix | Properties | Condition |
|--------|------------|-----------|
| Diagonal | Easy | κ = max/min diagonal |
| Random well-cond | General | κ ~ 10¹-10² |
| Hilbert | Ill-conditioned | κ ~ 10^n |
| Vandermonde | Ill-conditioned | κ grows fast |

## Interpreting Results

### Error Order Table

```
                sin    cos    exp    poly
CentralDiff      -8     -8     -8     -12
NDer4           -12    -12    -12    -16
NDer6           -14    -14    -14    -16
```

Reading:
- Row = algorithm
- Column = test function
- Value = log₁₀(error), more negative = better
- -16 is machine precision (best possible)

### Red Flags

1. **Error > 10⁻⁸** for smooth functions → Implementation bug
2. **Error varies wildly** across similar functions → Numerical instability
3. **Error doesn't improve** with smaller h → Wrong convergence order
4. **Error gets worse** with smaller h → Roundoff domination

### Expected Behaviors

| Method | Error vs h | Error vs n |
|--------|------------|------------|
| Central diff (NDer2) | O(h²) | N/A |
| NDer4/6/8 | O(h^(4/6/8)) | N/A |
| Simpson | N/A | O(n⁻⁴) |
| Gauss-Legendre | N/A | Exponential |
| RK4 | O(h⁴) | N/A |

## Timing Considerations

### Warm-up
```cpp
// Discard first call (cache warming)
compute(params);

// Now time
auto start = high_resolution_clock::now();
double result = compute(params);
auto end = high_resolution_clock::now();
```

### Multiple Runs
```cpp
const int RUNS = 100;
auto start = high_resolution_clock::now();
for (int i = 0; i < RUNS; i++) {
    result = compute(params);
}
auto end = high_resolution_clock::now();
double avg_time = duration<double, milli>(end - start).count() / RUNS;
```

### Avoid Optimization Artifacts
```cpp
volatile double result;  // Prevent dead code elimination
result = compute(params);
```

## CSV Export Format

```csv
Algorithm,Function,Exact,Computed,AbsoluteError,RelativeError,TimeMs,ErrorOrder
CentralDiff,sin(1),0.5403023058681398,0.5403023058681390,8.881784197001252e-16,1.643e-15,0.001,-15
NDer4,sin(1),0.5403023058681398,0.5403023058681398,0.000000000000000e+00,0.000e+00,0.003,-16
```

Use for:
- Plotting with external tools
- Regression testing
- Statistical analysis
- Documentation tables

## Adding New Tests

### Step 1: Plan Test Cases

- What algorithms to test?
- What functions/problems have exact solutions?
- What parameter ranges?
- What edge cases?

### Step 2: Implement Test Function

```cpp
void Test_Precision_NewFeature(MML::Tests::TestRegister& tests) {
    auto& suite = MML::Tests::GetPrecisionSuite();
    suite.clearResults();
    
    // Implement tests...
    
    suite.printSummary();
    suite.exportCSV("results/precision_newfeature.csv");
}
```

### Step 3: Register Test

In `main.cpp`:
```cpp
extern void Test_Precision_NewFeature(MML::Tests::TestRegister&);
// ...
Test_Precision_NewFeature(tests);
```

### Step 4: Update CMakeLists.txt

```cmake
set(TEST_SOURCES
    # ...
    src/testing_precision/test_precision_newfeature.cpp
)
```

### Step 5: Run and Validate

```bash
cmake --build build
cd build && ctest -V
```

### Step 6: Document Results

Create `docs/testing_precision/NEWFEATURE_ANALYSIS.md` with findings.

## Best Practices Checklist

- [ ] Use exact analytical solutions when possible
- [ ] Test multiple algorithms on same problem
- [ ] Vary parameters (h, n, tolerance) systematically
- [ ] Include ill-conditioned cases
- [ ] Verify convergence order
- [ ] Time performance (but accuracy first!)
- [ ] Export results for reproducibility
- [ ] Document unexpected behaviors
- [ ] Compare against theoretical bounds
- [ ] Test edge cases and boundary conditions

## Common Pitfalls

### 1. Testing at Only One Point
```cpp
// Bad: Only tests x = 1
suite.addResult("Method", "f(1)", exact, computed, time);

// Good: Test multiple points
for (double x : {0.1, 0.5, 1.0, 2.0, 10.0}) {
    // Test at each x
}
```

### 2. Ignoring Condition Number
```cpp
// Bad: "This solver failed!"
// Good: "This solver achieves κ·ε accuracy as expected"
```

### 3. Chasing Machine Precision
Not all problems can achieve 10⁻¹⁶. Know your problem's limitations.

### 4. Forgetting Roundoff Regime
```cpp
// For derivatives, optimal h ≈ ε^(1/3) ≈ 10⁻⁵
// Smaller h → roundoff dominates
```

### 5. Timing Noise
Single-run timings are noisy. Use averages for performance comparisons.

## Integration with CI/CD

### Regression Testing

```bash
# Run precision tests
./build/MML_Tests --precision

# Compare CSV outputs
python compare_precision.py results/old.csv results/new.csv

# Fail if accuracy degraded
```

### Automated Validation

```cpp
// Assert minimum accuracy
double error = std::abs(computed - exact);
REQUIRE(error < 1e-10);  // Catch2 assertion
```

## Conclusion

Precision testing is not optional for numerical libraries. The framework provides:

1. **Systematic methodology** for validating algorithms
2. **Quantitative metrics** for comparing approaches
3. **Documentation** of expected accuracy
4. **Regression detection** during development
5. **User guidance** for algorithm selection

Every numerical algorithm should have corresponding precision tests that verify it achieves its theoretical accuracy limits.

---

*Framework header: `include/core/PrecisionTestFramework.h`*
