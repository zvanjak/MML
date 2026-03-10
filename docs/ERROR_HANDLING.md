# MML Error Handling Policy

This document defines when and how MML code reports errors. It establishes the
two-tier philosophy used throughout the library and serves as the authoritative
reference for contributors.

> **See also:** [CODING_STYLE.md](CODING_STYLE.md) · [API_PATTERNS.md](API_PATTERNS.md) · [API_AUDIT_REPORT.md](API_AUDIT_REPORT.md)

---

## Table of Contents

1. [Two-Tier Philosophy](#two-tier-philosophy)
2. [Decision Flowchart](#decision-flowchart)
3. [Layer Rules](#layer-rules)
4. [Exception Hierarchy](#exception-hierarchy)
5. [Result Struct Anatomy](#result-struct-anatomy)
6. [Writing Good Error Messages](#writing-good-error-messages)
7. [Caller Patterns](#caller-patterns)
8. [Known Inconsistencies](#known-inconsistencies)
9. [Migration Guide: Throw → Result](#migration-guide-throw--result)

---

## Two-Tier Philosophy

MML uses **two distinct mechanisms** for reporting errors, chosen by the
*nature* of the failure:

| Mechanism | When to Use | Example |
|-----------|------------|---------|
| **Exceptions** | Programmer errors — precondition violations that indicate a bug in calling code | Dimension mismatch, index out-of-bounds, null function, negative tolerance |
| **Result structs** | Expected algorithmic failures — situations that can arise from valid inputs | No convergence, singular matrix detected at runtime, step-size underflow |

### Why Two Mechanisms?

- **Exceptions** propagate up the call stack automatically. They are appropriate
  for conditions the caller *must never cause* — recovering from them is usually
  impossible without code changes.

- **Result structs** keep control flow local. They are appropriate for conditions
  the caller *should anticipate and handle* — e.g., an iterative solver that
  fails to converge with the given tolerance should report that fact, not crash.

### The Golden Rules

1. **Never throw for expected algorithmic failure.** If a solver can reasonably
   fail to converge, return a Result with `converged = false`.

2. **Never return a Result for programmer error.** If the caller passes a 3×4
   matrix to a function requiring square matrices, throw `MatrixDimensionError`.

3. **Never swallow exceptions in algorithm code.** If a lower-level operation
   throws (e.g., singular matrix during LU decomposition inside a Newton step),
   catch it and convert to a failed Result — do not let it propagate as if it
   were a precondition violation of the algorithm itself.

4. **Never throw *and* return an error for the same failure.** Pick one. Algorithms
   that populate a Result and then `throw` the error message are an anti-pattern.

---

## Decision Flowchart

```
Is the failure caused by invalid arguments that the caller should have validated?
  YES → throw MML exception (precondition violation)
  NO  ↓

Is the failure an expected outcome for some valid inputs?
(e.g., no convergence, singular matrix encountered during iteration)
  YES → return Result struct with converged=false / status != Success
  NO  ↓

Is this a system-level failure? (out of memory, file I/O)
  YES → throw std::runtime_error or MML::FileIOError
  NO  → consult the team
```

---

## Layer Rules

Each MML layer has a primary error-handling mechanism:

### Base Layer (`mml/base/`)

**Primary:** Exceptions only.

Vector, Matrix, Polynom, Tensor, and other fundamental types validate
preconditions and throw on violation. These types never return Result structs.

```cpp
// Matrix.h — precondition check
if (i < 0 || i >= rows())
    throw MatrixAccessBoundsError("Matrix::operator()", i, j, rows(), cols());
```

**Rationale:** Base types are building blocks with strict contracts. Any misuse
is a programmer error.

### Core Layer (`mml/core/`)

**Primary:** Exceptions for direct operations, Result structs for iterative solvers.

- **Direct solvers** (LU, QR, SVD, Cholesky) validate inputs with exceptions,
  and may throw `SingularMatrixError` when the *input itself* is invalid (e.g.,
  the caller was supposed to check for singularity before calling).
- **Iterative solvers** (Jacobi iteration, Gauss-Seidel, SOR) return
  `IterativeSolverResult` with convergence status.
- **Integration** returns `IntegrationResult` structs.
- **Derivation** throws for invalid inputs (e.g., step size ≤ 0).

### Algorithms Layer (`mml/algorithms/`)

**Primary:** Result structs for all iterative/optimization/root-finding algorithms.

Every algorithm that iterates toward a solution **must** return a Result struct.
The Result should include at minimum: `converged`, `iterations_used`,
`achieved_tolerance`, and `error_message`.

Exceptions are still appropriate for *precondition validation at function entry*:

```cpp
RootFindingResult FindRootBrent(const IRealFunction& f, Real a, Real b,
                                 const RootFindingConfig& config) {
    // Precondition check — throw
    if (config.tolerance <= 0)
        throw ArgumentError("FindRootBrent: tolerance must be positive");

    // Algorithm may fail to converge — return Result
    RootFindingResult result;
    result.algorithm_name = "Brent";
    // ... iteration ...
    if (iter >= config.max_iterations) {
        result.converged = false;
        result.status = AlgorithmStatus::MaxIterationsExceeded;
        result.error_message = "Failed to converge after " +
                               std::to_string(iter) + " iterations";
        return result;  // DO NOT throw here
    }
}
```

**Reference implementation:** [`mml/algorithms/RootFinding.h`](../mml/algorithms/RootFinding.h)

### Systems Layer (`mml/systems/`)

**Primary:** Follow the same rules as Algorithms — Result structs for iterative
operations, exceptions for precondition violations.

### Tools Layer (`mml/tools/`)

**Primary:** Result structs with success/failure factories.

Tools (serializer, data loader, task runner) use their own Result types with
`bool success`, `std::string errorMessage`, and optional `warningMessages`.
They catch internal exceptions and translate them into Result fields.

```cpp
// Serializer pattern
auto result = SerializeResult::Failure("Failed to open file: " + path);
auto result = SerializeResult::Success("Wrote 42 items");
```

### Packages Layer (`mml_packages/`)

**Primary:** Follow the Algorithms layer rules — Result structs for iterative
operations.

---

## Exception Hierarchy

All MML exceptions are defined in [`mml/MMLExceptions.h`](../mml/MMLExceptions.h).
They inherit from standard library exception types:

### Generic Errors

| Exception | Base Class | When to Throw |
|-----------|-----------|---------------|
| `ArgumentError` | `std::invalid_argument` | Invalid argument value |
| `DomainError` | `std::domain_error` | Mathematical domain violation |
| `DivisionByZeroError` | `std::domain_error` | Division by zero |
| `NotImplementedError` | `std::logic_error` | Abstract method / not-yet-implemented |
| `IndexError` | `std::out_of_range` | Generic index out of range |

### Data Structure Errors

| Exception | Base Class | When to Throw |
|-----------|-----------|---------------|
| `VectorInitializationError` | `std::invalid_argument` | Invalid Vector construction |
| `VectorDimensionError` | `std::invalid_argument` | Vector size mismatch |
| `VectorAccessBoundsError` | `std::out_of_range` | Vector index out of range |
| `MatrixAllocationError` | `std::out_of_range` | Invalid Matrix allocation |
| `MatrixAccessBoundsError` | `std::out_of_range` | Matrix index out of range |
| `MatrixDimensionError` | `std::invalid_argument` | Matrix dimension mismatch |
| `SingularMatrixError` | `std::domain_error` | Singular or near-singular matrix |
| `MatrixNumericalError` | `std::runtime_error` | Numerical failure in matrix operation |

### Tensor Errors

| Exception | Base Class | When to Throw |
|-----------|-----------|---------------|
| `TensorCovarContravarNumError` | `std::invalid_argument` | Wrong index count |
| `TensorCovarContravarArithmeticError` | `std::invalid_argument` | Invalid tensor arithmetic |
| `TensorIndexError` | `std::invalid_argument` | Invalid tensor index |

### Algorithm / Numerical Errors

| Exception | Base Class | When to Throw |
|-----------|-----------|---------------|
| `NumericalMethodError` | `std::runtime_error` | Generic numerical failure |
| `ConvergenceError` | `std::runtime_error` | Convergence failure (**prefer Result struct**) |
| `RootFindingError` | `std::runtime_error` | Root finding failure (**prefer Result struct**) |
| `ODESolverError` | `std::runtime_error` | ODE solver failure (**prefer Result struct**) |
| `StatisticsError` | `std::runtime_error` | Statistics computation failure |
| `FourierError` | `std::invalid_argument` | Invalid Fourier input (non-power-of-2) |

### Domain-Specific Errors

| Exception | Base Class | When to Throw |
|-----------|-----------|---------------|
| `IntegrationTooManySteps` | `std::domain_error` | Integration exceeded step limit |
| `RealFuncInterpInitError` | `std::domain_error` | Interpolation setup failure |
| `RealFuncInterpRuntimeError` | `std::runtime_error` | Interpolation evaluation failure |
| `GeometryError` | `std::domain_error` | Invalid geometric operation |
| `QuaternionError` | `std::domain_error` | Invalid quaternion operation |
| `CurveFittingError` | `std::invalid_argument` | Invalid curve fitting input |
| `DataError` | `std::runtime_error` | Data validation failure |
| `NumericInputError` | `std::domain_error` | Invalid numeric input |

### I/O and Tools Errors

| Exception | Base Class | When to Throw |
|-----------|-----------|---------------|
| `FileIOError` | `std::runtime_error` | File read/write failure |
| `VisualizerError` | `std::runtime_error` | Visualization failure |

### Important Note on Algorithm Exceptions

`ConvergenceError`, `RootFindingError`, and `ODESolverError` exist for
**legacy compatibility**. New algorithm code should **not** throw these —
use Result structs instead. These exceptions remain available for callers who
prefer exceptions via opt-in wrapper functions.

---

## Result Struct Anatomy

The canonical base pattern is defined in
[`mml/core/AlgorithmTypes.h`](../mml/core/AlgorithmTypes.h):

### Required Fields

Every algorithm Result struct **must** include these fields:

```cpp
struct MyAlgorithmResult {
    // === Convergence status ===
    bool converged = false;                           // Did it succeed?
    int iterations_used = 0;                          // How many iterations?
    Real achieved_tolerance = 0.0;                    // What precision was achieved?
    AlgorithmStatus status = AlgorithmStatus::Success; // Structured status code
    std::string error_message;                        // Empty on success, descriptive on failure

    // === Algorithm output ===
    // ... algorithm-specific fields (root, eigenvalues, solution, etc.)
};
```

### Recommended Fields

```cpp
    std::string algorithm_name;     // "Brent", "Jacobi", "RKF45"
    double elapsed_time_ms = 0.0;   // Wall-clock time
    int function_evaluations = 0;   // Number of f(x) calls
```

### AlgorithmStatus Codes

| Code | Meaning |
|------|---------|
| `Success` | Algorithm converged successfully |
| `MaxIterationsExceeded` | Hit iteration limit without converging |
| `NumericalInstability` | NaN, Inf, or ill-conditioning detected |
| `SingularMatrix` | Singular matrix encountered during iteration |
| `InvalidInput` | Input validation failed |
| `Stalled` | No progress between iterations |
| `ToleranceUnachievable` | Requested tolerance cannot be reached |
| `AlgorithmSpecificFailure` | See `error_message` for details |

### Reference Implementations

- **Root Finding:** `RootFindingResult` in [`mml/algorithms/RootFinding.h`](../mml/algorithms/RootFinding.h)
- **Eigen Solvers:** `EigenSolverResult` in [`mml/algorithms/EigenSystemSolvers.h`](../mml/algorithms/EigenSystemSolvers.h)
- **ODE Solvers:** `ODEIntegrationResult` in [`mml/algorithms/ODESolvers/ODESystemSteppers.h`](../mml/algorithms/ODESolvers/ODESystemSteppers.h)
- **Optimization:** `OptimizationResult` in [`mml/algorithms/Optimization.h`](../mml/algorithms/Optimization.h)
- **Tools:** `SerializeResult` in [`mml/tools/serializer/SerializerBase.h`](../mml/tools/serializer/SerializerBase.h)

---

## Writing Good Error Messages

### Exception Messages

Include **what failed**, **why**, and **relevant values**:

```cpp
// ✅ Good — specific, actionable
throw MatrixDimensionError(
    "Matrix multiply: A(3×4) × B(5×2) — inner dimensions don't match",
    3, 4, 5, 2);

// ❌ Bad — vague
throw std::runtime_error("dimension error");
```

### Result Error Messages

Include **what the algorithm tried**, **how far it got**, and **what stopped it**:

```cpp
// ✅ Good — diagnostic-rich
result.error_message =
    "Brent: failed to converge after 100 iterations "
    "(residual: 1.5e-4, bracket: [1.41421, 1.41422])";

// ❌ Bad — useless
result.error_message = "did not converge";
```

---

## Caller Patterns

### Checking a Result

```cpp
auto result = RootFinding::FindRootBrent(f, 0.0, 2.0, config);

if (result.converged) {
    // Safe to use result.root
    std::cout << "Root: " << result.root << "\n";
} else {
    // Handle failure — log, retry with different parameters, fall back
    std::cerr << "Warning: " << result.error_message << "\n";
}
```

### Catching Exceptions

```cpp
try {
    auto LU = A.LUDecompose();
    auto x = LU.Solve(b);
} catch (const MatrixDimensionError& e) {
    // Programmer error — fix the calling code
    std::cerr << "Bug: " << e.what() << "\n";
} catch (const SingularMatrixError& e) {
    // Possibly expected — handle gracefully
    std::cerr << "Matrix is singular: " << e.what() << "\n";
}
```

### Converting Between Tiers

If higher-level code wants exceptions from an algorithm that returns Result:

```cpp
auto result = solver.Solve(system, config);
if (!result.converged) {
    throw ConvergenceError("Solver failed: " + result.error_message);
}
```

If algorithm internals encounter an exception that should become a Result:

```cpp
try {
    // LU solve inside Newton iteration
    auto delta = jacobian.Solve(residual);
} catch (const SingularMatrixError&) {
    result.converged = false;
    result.status = AlgorithmStatus::SingularMatrix;
    result.error_message = "Jacobian became singular at iteration " +
                           std::to_string(iter);
    return result;  // Convert exception → Result
}
```

---

## Known Inconsistencies

The following areas do not yet fully conform to this policy. They are tracked
for future cleanup but do not affect correctness.

### Algorithms That Throw for Expected Failures

These functions throw exceptions where they should return a failed Result:

| File | Function | Throws | Should Return |
|------|----------|--------|---------------|
| `RootFindingPolynoms.h` | `LaguerreRoot` | `ConvergenceError` | `RootFindingResult` with `converged=false` |
| `RootFindingBracketing.h` | Legacy wrappers | Throws **after** populating Result | Just return the Result |
| `Optimization.h` | `GoldenSectionMin` | `ConvergenceError` on bracket failure | `OptimizationResult` with `converged=false` |
| `Optimization.h` | `BrentMin` | `ConvergenceError` on max iterations | `OptimizationResult` with `converged=false` |
| `ODEStiffSolvers.h` | `BackwardEuler`, `BDF2` | `ODESolverError` on Newton failure | Return via `ODEIntegrationResult` |
| `ODEStiffSolvers.h` | `Rosenbrock23` | `ODESolverError` on linear solve failure | Return via `ODEIntegrationResult` |
| `ODEAdaptiveIntegrator.h` | Adaptive stepper | `ODESolverError` on step-size underflow | Return via `ODEIntegrationResult` |

### Raw `std::` Exceptions Instead of MML Types

| File | Uses | Should Use |
|------|------|------------|
| `ConvexHull.h` | `std::runtime_error`, `std::invalid_argument` | `GeometryError`, `ArgumentError` |
| `GraphSpectral.h` | `std::runtime_error` | `ArgumentError` or `NumericalMethodError` |

### Result Structs Missing Standard Fields

Some older Result structs predate the `IterativeResultBase` pattern and lack
`AlgorithmStatus`, `achieved_tolerance`, or `error_message`. These are functional
but provide less diagnostic information than newer Result types.

---

## Migration Guide: Throw → Result

When converting an algorithm from throwing to returning a Result:

### Step 1: Identify the Failure Points

Find all `throw` statements that represent **expected algorithmic failures**
(not precondition violations):

```cpp
// BEFORE: throws on convergence failure
if (iter >= maxIter)
    throw ConvergenceError("Failed to converge");
```

### Step 2: Replace with Result Population

```cpp
// AFTER: returns failed Result
if (iter >= maxIter) {
    result.converged = false;
    result.status = AlgorithmStatus::MaxIterationsExceeded;
    result.iterations_used = iter;
    result.error_message = "Failed to converge after " +
                           std::to_string(iter) + " iterations"
                           " (residual: " + std::to_string(residual) + ")";
    return result;
}
```

### Step 3: Ensure Success Path Populates Result

```cpp
// On success
result.converged = true;
result.status = AlgorithmStatus::Success;
result.iterations_used = iter;
result.achieved_tolerance = residual;
return result;
```

### Step 4: Update Tests

Tests that used `REQUIRE_THROWS` for convergence failures should switch to
checking the Result:

```cpp
// BEFORE
REQUIRE_THROWS_AS(solver.Solve(badInput), ConvergenceError);

// AFTER
auto result = solver.Solve(badInput, config);
REQUIRE_FALSE(result.converged);
REQUIRE(result.status == AlgorithmStatus::MaxIterationsExceeded);
```

### Step 5: Keep Precondition Throws

Do **not** convert precondition checks to Result. These should remain as throws:

```cpp
// KEEP THIS — it's a programmer error
if (A.rows() != A.cols())
    throw MatrixDimensionError("Requires square matrix", ...);
```

---

*Document Version: 1.0*  
*Last Updated: 2026-02-07*  
*Task: MinimalMathLibrary-2l5d.25 (Document error handling policy)*
