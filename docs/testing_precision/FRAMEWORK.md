# PrecisionTestFramework - Design and API

## Overview

The `PrecisionTestFramework.h` provides a unified infrastructure for precision testing numerical algorithms. It enables standardized result collection, statistical analysis, and multiple output formats.

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                    PrecisionTestSuite                           │
│  ┌───────────────────────────────────────────────────────────┐  │
│  │                PrecisionTestResult                        │  │
│  │  • algorithm_name    • exact_value                        │  │
│  │  • test_function     • computed_value                     │  │
│  │  • parameters        • absolute_error                     │  │
│  │  • iterations        • relative_error                     │  │
│  │  • time_ms           • converged                          │  │
│  │  • notes             • getErrorOrder()                    │  │
│  └───────────────────────────────────────────────────────────┘  │
│                                                                 │
│  ┌───────────────────────────────────────────────────────────┐  │
│  │               PrecisionStatistics                         │  │
│  │  • min/max/avg absolute error                             │  │
│  │  • min/max/avg relative error                             │  │
│  │  • best/worst error order                                 │  │
│  │  • convergence rate                                       │  │
│  └───────────────────────────────────────────────────────────┘  │
│                                                                 │
│  Methods:                                                       │
│  • addResult()          • printDetailedTable()                  │
│  • calculateStatistics()• printErrorOrderMatrix()               │
│  • printSummary()       • exportCSV()                           │
└─────────────────────────────────────────────────────────────────┘
```

## Core Classes

### PrecisionTestResult

Captures a single precision test result with all relevant metrics.

```cpp
struct PrecisionTestResult {
    std::string algorithm_name;     // Name of the algorithm
    std::string test_function;      // Name of the test function
    std::string parameters;         // Algorithm parameters (e.g., "h=1e-6")
    
    double exact_value;             // Known exact/analytical value
    double computed_value;          // Value computed by the algorithm
    double absolute_error;          // |computed - exact|
    double relative_error;          // |computed - exact| / |exact|
    
    int iterations;                 // Number of iterations (if applicable)
    double time_ms;                 // Computation time in milliseconds
    
    bool converged;                 // Whether the algorithm converged
    std::string notes;              // Additional notes or warnings
    
    int getErrorOrder() const;      // Returns floor(log10(absolute_error))
};
```

### PrecisionStatistics

Aggregated statistics for a collection of test results.

```cpp
struct PrecisionStatistics {
    double min_abs_error, max_abs_error, avg_abs_error, stddev_abs_error;
    double min_rel_error, max_rel_error, avg_rel_error, stddev_rel_error;
    
    int total_tests;
    int converged_tests;
    double avg_time_ms;
    
    int best_error_order;           // Best (most negative) error order
    int worst_error_order;          // Worst (least negative) error order
};
```

### PrecisionTestSuite

Container and analyzer for test results.

```cpp
class PrecisionTestSuite {
public:
    PrecisionTestSuite(const std::string& name, const std::string& desc = "");
    
    // Adding results
    void addResult(const PrecisionTestResult& result);
    void addResult(const std::string& algo, const std::string& func,
                   double exact, double computed, double time_ms = 0.0);
    
    // Statistics
    PrecisionStatistics calculateStatistics() const;
    PrecisionStatistics calculateStatistics(const std::string& algorithm) const;
    
    // Output
    void printHeader(std::ostream& os = std::cout) const;
    void printDetailedTable(std::ostream& os = std::cout) const;
    void printErrorOrderMatrix(std::ostream& os = std::cout) const;
    void printSummary(std::ostream& os = std::cout) const;
    
    // Export
    void exportCSV(const std::string& filename) const;
    void exportMarkdown(const std::string& filename) const;
};
```

## Usage Patterns

### Basic Test Pattern

```cpp
void Test_MyAlgorithm(PrecisionTestSuite& suite) {
    // Known exact value
    double exact = 2.718281828459045;  // e
    
    // Compute using algorithm
    double computed = myAlgorithm();
    
    // Record result
    suite.addResult("MyAlgorithm", "exp(1)", exact, computed);
}
```

### Timing Pattern

```cpp
void Test_WithTiming(PrecisionTestSuite& suite) {
    double exact = 3.141592653589793;
    
    auto start = std::chrono::high_resolution_clock::now();
    double computed = computePi();
    auto end = std::chrono::high_resolution_clock::now();
    
    double time_ms = std::chrono::duration<double, std::milli>(end - start).count();
    
    suite.addResult("ComputePi", "pi", exact, computed, time_ms);
}
```

### Parameter Sweep Pattern

```cpp
void Test_ParameterSweep(PrecisionTestSuite& suite) {
    double exact = 2.0;  // f'(x) = 2x at x=1
    
    for (double h : {1e-2, 1e-4, 1e-6, 1e-8, 1e-10}) {
        double computed = centralDifference(f, 1.0, h);
        
        PrecisionTestResult r("CentralDiff", "x^2", exact, computed);
        r.parameters = "h=" + std::to_string(h);
        suite.addResult(r);
    }
}
```

### Convergence Order Pattern

```cpp
void Test_ConvergenceOrder(PrecisionTestSuite& suite) {
    std::vector<double> errors;
    std::vector<double> step_sizes = {0.1, 0.05, 0.025, 0.0125};
    
    for (double h : step_sizes) {
        double computed = integrate(f, 0, 1, h);
        double error = std::abs(computed - exact);
        errors.push_back(error);
        
        suite.addResult("Trapezoid", "sin(x)", exact, computed);
    }
    
    // Calculate convergence order
    for (size_t i = 1; i < errors.size(); i++) {
        double order = std::log(errors[i-1]/errors[i]) / std::log(2.0);
        std::cout << "Order: " << order << std::endl;  // Should be ~2 for trapezoid
    }
}
```

## Output Formats

### Error Order Matrix

The `printErrorOrderMatrix()` method produces a compact view showing the order of magnitude of errors for each algorithm/function combination:

```
ERROR ORDER MATRIX (log10 of absolute error):
--------------------------------------------------------------------------------
Algorithm      f(x)=x²   f(x)=sin  f(x)=exp  f(x)=log  AVG
--------------------------------------------------------------------------------
NDer1           -8        -8        -8        -7      -7.8
NDer2          -12       -12       -12       -11     -11.8
NDer6          -14       -14       -14       -13     -13.8
--------------------------------------------------------------------------------
```

### CSV Export

The `exportCSV()` method produces a comma-separated file for further analysis:

```csv
Algorithm,Function,Parameters,Exact,Computed,AbsError,RelError,Order,Time_ms,Converged
NDer2,x^2,h=1e-6,2.000000,2.000000,1.23e-12,6.15e-13,-12,0.001,true
NDer6,x^2,h=1e-4,2.000000,2.000000,4.56e-14,2.28e-14,-14,0.002,true
```

## Best Practices

### 1. Test Function Selection

Choose test functions that:
- Have known exact solutions
- Cover different difficulty levels (smooth, oscillatory, singular)
- Exercise edge cases (zeros, infinities, discontinuities)

### 2. Error Interpretation

```
Order -16: Machine precision (excellent)
Order -14 to -15: Near-optimal
Order -10 to -13: Good for most applications
Order -6 to -9: Acceptable for visualization/quick estimates
Order > -6: May indicate algorithm issues or ill-conditioning
```

### 3. Convergence Verification

Always verify that iterative methods have converged:

```cpp
PrecisionTestResult r("Newton", "x^3-2", exact, computed);
r.converged = (iterations < max_iterations);
r.iterations = iterations;
if (!r.converged) {
    r.notes = "Did not converge within " + std::to_string(max_iterations) + " iterations";
}
suite.addResult(r);
```

### 4. Comparative Analysis

Use the framework to compare algorithms:

```cpp
// Run same test with multiple algorithms
suite.addResult("Bisection", "sqrt(2)", exact, bisection(f, 1, 2));
suite.addResult("Newton", "sqrt(2)", exact, newton(f, 1.5));
suite.addResult("Brent", "sqrt(2)", exact, brent(f, 1, 2));

// Get per-algorithm statistics
auto bisection_stats = suite.calculateStatistics("Bisection");
auto newton_stats = suite.calculateStatistics("Newton");
auto brent_stats = suite.calculateStatistics("Brent");
```

## Extension Points

The framework is designed for extension:

1. **Custom output formats** - Implement new export methods
2. **Additional statistics** - Extend PrecisionStatistics
3. **Visualization** - Generate data for plotting convergence curves
4. **Regression testing** - Compare results against baseline

---

*See individual algorithm analysis documents for detailed findings.*
