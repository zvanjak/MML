# MML API Patterns Reference

This document defines the standard patterns for all MML algorithm APIs. Follow these patterns when creating new algorithms or updating existing ones.

## Table of Contents

1. [Overview](#overview)
2. [Config Objects](#config-objects)
3. [Result Structs](#result-structs)
4. [Function Signatures](#function-signatures)
5. [Error Handling](#error-handling)
6. [Backward Compatibility](#backward-compatibility)
7. [Examples](#examples)

---

## Overview

MML uses a consistent **Config + Result** pattern across all iterative algorithms:

```cpp
// User configures the algorithm
AlgorithmConfig config;
config.tolerance = 1e-12;
config.max_iterations = 200;

// Algorithm returns a rich result
AlgorithmResult result = SomeAlgorithm(input, config);

// User checks convergence and accesses data
if (result.converged) {
    std::cout << "Solution: " << result.value << "\n";
    std::cout << "Achieved in " << result.iterations_used << " iterations\n";
}
```

### Benefits

- **Self-documenting**: Named fields instead of positional parameters
- **Extensible**: Add fields without breaking API
- **Diagnostic-rich**: Users get convergence info, not just raw values
- **Testable**: Easy to verify algorithm behavior
- **IDE-friendly**: Autocomplete shows all options

### Reference Implementation

**`mml/algorithms/RootFinding.h`** is the reference implementation. When in doubt, follow its patterns.

---

## Config Objects

### Standard Structure

```cpp
/// Configuration for [Algorithm Family] algorithms.
/// 
/// Provides user control over convergence criteria, iteration limits,
/// and diagnostic output. All fields have sensible defaults.
struct AlgorithmConfig {
    /// Convergence tolerance (default: 1e-10)
    Real tolerance = 1e-10;
    
    /// Maximum iterations (default: 100)
    int max_iterations = 100;
    
    /// Enable verbose output (default: false)
    bool verbose = false;
    
    // Algorithm-specific fields below...
};
```

### Required Fields

All Config objects should include:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `tolerance` | `Real` | `1e-10` | Convergence tolerance |
| `max_iterations` | `int` | `100` | Maximum iteration count |
| `verbose` | `bool` | `false` | Enable debug output |

### Optional Fields (by algorithm family)

**Eigen Solvers:**
- `sort_eigenvalues` (bool) - Sort eigenvalues by magnitude
- `compute_eigenvectors` (bool) - Compute eigenvectors (some methods can skip)

**ODE Integrators:**
- `initial_step_size` (Real) - Starting step size
- `min_step_size` (Real) - Minimum allowed step
- `max_step_size` (Real) - Maximum allowed step
- `abs_tolerance` (Real) - Absolute error tolerance
- `rel_tolerance` (Real) - Relative error tolerance
- `dense_output` (bool) - Store intermediate points

**Optimization:**
- `gradient_tolerance` (Real) - Tolerance for gradient norm
- `line_search_config` (nested) - Line search parameters

**Iterative Linear Solvers:**
- `omega` (Real) - Relaxation factor for SOR
- `preconditioner_type` (enum) - Preconditioner selection

### Naming Convention

- Config struct: `{AlgorithmFamily}Config` (e.g., `EigenSolverConfig`, `ODEIntegratorConfig`)
- All field names: `snake_case`
- Boolean fields: Positive phrasing (e.g., `verbose` not `quiet`)

---

## Result Structs

### Standard Structure

```cpp
/// Result of [Algorithm] operation.
/// 
/// Contains the solution along with diagnostic information about
/// the convergence process. Always check `converged` before using the solution.
struct AlgorithmResult {
    // === Primary Output ===
    /// The computed solution
    SolutionType solution;
    
    // === Convergence Status ===
    /// True if algorithm converged within tolerance
    bool converged = false;
    
    /// Number of iterations used
    int iterations_used = 0;
    
    /// Achieved tolerance (may be better than requested)
    Real achieved_tolerance = 0.0;
    
    // === Diagnostics ===
    /// Name of the algorithm that produced this result
    std::string algorithm_name;
    
    /// Error message if not converged (empty on success)
    std::string error_message;
    
    // === Performance (optional) ===
    /// Elapsed time in milliseconds
    double elapsed_time_ms = 0.0;
    
    /// Number of function evaluations
    int function_evaluations = 0;
};
```

### Required Fields

All Result structs MUST include:

| Field | Type | Description |
|-------|------|-------------|
| `converged` | `bool` | Did algorithm succeed? |
| `iterations_used` | `int` | Actual iteration count |
| `achieved_tolerance` | `Real` | Actual achieved accuracy |
| `error_message` | `std::string` | Error description (empty on success) |

### Recommended Fields

| Field | Type | Description |
|-------|------|-------------|
| `algorithm_name` | `std::string` | Algorithm identifier (for logging) |
| `elapsed_time_ms` | `double` | Wall-clock time in ms |
| `function_evaluations` | `int` | Number of f(x) calls |

### Algorithm-Specific Fields

**Root Finding:**
- `root` (Real) - The computed root
- `function_value` (Real) - f(root), should be ~0

**Eigen Solvers:**
- `eigenvalues` (Vector) - Computed eigenvalues
- `eigenvectors` (Matrix) - Computed eigenvectors
- `max_residual` (Real) - Maximum |Av - λv|

**ODE Integrators:**
- `final_state` (Vector) - State at final time
- `final_time` (Real) - Actual end time
- `total_steps` (int) - Number of steps taken
- `rejected_steps` (int) - Steps rejected for error control

**Optimization:**
- `minimum_value` (Real) - The minimum f(x) found
- `minimizer` (Vector) - x that minimizes f
- `gradient_norm` (Real) - |∇f| at solution

### Naming Convention

- Result struct: `{AlgorithmFamily}Result` (e.g., `EigenSolverResult`, `ODEIntegrationResult`)
- All field names: `snake_case`
- Primary output field: Named for what it contains (`root`, `solution`, `eigenvalues`)

---

## Function Signatures

### Standard Pattern

```cpp
// Config-based (preferred)
AlgorithmResult Algorithm(const InputType& input, const AlgorithmConfig& config);

// Simple (for backward compatibility, delegates to config-based)
AlgorithmResult Algorithm(const InputType& input, Real tol = 1e-10, int max_iter = 100);
```

### Parameter Passing Rules

| Parameter Type | Pass By | Example |
|---------------|---------|---------|
| Input matrices/vectors | `const Type&` | `const Matrix<Real>& A` |
| Input scalars | Value | `Real tolerance` |
| Config objects | `const Config&` | `const EigenSolverConfig& config` |
| Output (return) | By value | `EigenSolverResult` |
| Optional outputs | Pointer with nullptr | `Matrix<Real>* eigenvectors = nullptr` |

### Example Signatures

```cpp
// Root finding
RootFindingResult FindRootBrent(const IRealFunction& func, 
                                 Real x1, Real x2,
                                 const RootFindingConfig& config);

// Eigen solver
EigenSolverResult Jacobi_Solve(const Matrix<Real>& A,
                                const EigenSolverConfig& config);

// ODE integrator
ODEIntegrationResult Integrate(const IODESystem<Real>& system,
                                const Vector<Real>& y0,
                                Real t_start, Real t_end,
                                const ODEIntegratorConfig& config);
```

---

## Error Handling

> **Full policy:** See [ERROR_HANDLING.md](ERROR_HANDLING.md) for the complete error
> handling policy, exception hierarchy, known inconsistencies, and the
> throw-to-Result migration guide.

### Philosophy

- **Expected failures** (no convergence, singular matrix): Return in Result struct
- **Programmer errors** (invalid dimensions, null pointer): Throw exception

### Result-Based Error Reporting

```cpp
struct AlgorithmResult {
    bool converged = false;
    std::string error_message;  // Non-empty if failed
    
    // Optional: structured error codes
    enum class Status {
        Success,
        MaxIterationsExceeded,
        NumericalInstability,
        SingularMatrix,
        InvalidInput
    };
    Status status = Status::Success;
};
```

### Error Messages

Error messages should be:
- **Specific**: "Failed to converge after 100 iterations (residual: 1.5e-4)"
- **Actionable**: Include relevant values for debugging
- **Empty on success**: `error_message.empty() == true` means success

### Exception Use Cases

Throw `MML::MathException` for:
- Matrix dimension mismatch
- Null function pointer
- Invalid configuration (negative tolerance)
- Out-of-bounds access

```cpp
if (A.RowNum() != A.ColNum()) {
    throw MathException("Jacobi_Solve requires square matrix, got " +
                       std::to_string(A.RowNum()) + "x" + std::to_string(A.ColNum()));
}
```

---

## Backward Compatibility

### Deprecated Conversion Operators

To maintain backward compatibility while encouraging proper usage:

```cpp
struct RootFindingResult {
    Real root;
    bool converged;
    
    /// Conversion for backward compatibility
    /// @deprecated Use .root instead - implicit conversion discards diagnostics
    [[deprecated("Use .root instead - implicit conversion discards convergence info")]]
    operator Real() const { return root; }
};
```

### Simple Overloads

Maintain simple function overloads that delegate to config-based versions:

```cpp
// Simple API (delegates to config-based)
inline RootFindingResult FindRootBrent(const IRealFunction& func, 
                                        Real x1, Real x2,
                                        Real tol = 1e-10) {
    RootFindingConfig config;
    config.tolerance = tol;
    return FindRootBrent(func, x1, x2, config);
}

// Full API with config
RootFindingResult FindRootBrent(const IRealFunction& func, 
                                 Real x1, Real x2,
                                 const RootFindingConfig& config);
```

---

## Examples

### Adding Config to an Existing Algorithm

**Before:**
```cpp
static void Jacobi_Solve(const Matrix<Real>& A,
                         Vector<Real>& eigenvalues,
                         Matrix<Real>& eigenvectors,
                         Real tolerance = 1e-10,
                         int max_iter = 100);
```

**After:**
```cpp
struct EigenSolverConfig {
    Real tolerance = 1e-10;
    int max_iterations = 100;
    bool sort_eigenvalues = true;
    bool compute_eigenvectors = true;
    bool verbose = false;
};

struct EigenSolverResult {
    Vector<Real> eigenvalues;
    Matrix<Real> eigenvectors;
    bool converged = false;
    int iterations_used = 0;
    Real max_residual = 0.0;
    std::string algorithm_name = "Jacobi";
    std::string error_message;
};

// New config-based API
static EigenSolverResult Jacobi_Solve(const Matrix<Real>& A,
                                       const EigenSolverConfig& config);

// Backward-compatible wrapper
static void Jacobi_Solve(const Matrix<Real>& A,
                         Vector<Real>& eigenvalues,
                         Matrix<Real>& eigenvectors,
                         Real tolerance = 1e-10,
                         int max_iter = 100) {
    EigenSolverConfig config;
    config.tolerance = tolerance;
    config.max_iterations = max_iter;
    auto result = Jacobi_Solve(A, config);
    eigenvalues = result.eigenvalues;
    eigenvectors = result.eigenvectors;
}
```

### Complete Working Example

```cpp
#include "mml/algorithms/RootFinding.h"

using namespace MML;

// Custom configuration
RootFinding::RootFindingConfig config;
config.tolerance = 1e-12;
config.max_iterations = 200;
config.verbose = true;

// Define function
RealFunction f([](Real x) { return x*x - 2; });

// Find root with diagnostics
auto result = RootFinding::FindRootBrent(f, 0, 2, config);

if (result.converged) {
    std::cout << "√2 = " << std::setprecision(15) << result.root << "\n";
    std::cout << "Iterations: " << result.iterations_used << "\n";
    std::cout << "Residual: " << result.function_value << "\n";
} else {
    std::cerr << "Error: " << result.error_message << "\n";
}
```

---

## Checklist for New Algorithms

When adding a new algorithm:

- [ ] Create `AlgorithmConfig` struct with standard fields
- [ ] Create `AlgorithmResult` struct with convergence info
- [ ] Implement config-based function
- [ ] Add backward-compatible simple overload
- [ ] Add `[[deprecated]]` conversion operator if result has single primary value
- [ ] Document all fields with Doxygen comments
- [ ] Add unit tests for config parameters
- [ ] Add test for convergence reporting
- [ ] Update this document if introducing new patterns

---

## See Also

- [RootFinding.h](../../mml/algorithms/RootFinding.h) - Reference implementation
- [AlgorithmTypes.h](../../mml/core/AlgorithmTypes.h) - Base types
- [PARAMETER_PASSING.md](PARAMETER_PASSING.md) - Parameter passing conventions

---

*Document Version: 1.0*  
*Last Updated: 2026-01-21*  
*Epic: MinimalMathLibrary-rynt (API Standardization Initiative)*
