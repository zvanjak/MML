# Thread-Safety Guide

This document describes the thread-safety guarantees for MinimalMathLibrary classes and functions.

## Quick Reference

| Category | Meaning | Example |
|----------|---------|---------|
| **Thread-Safe** | Safe for concurrent reads; writes need synchronization | `Vector::operator[] const` |
| **Reentrant** | Safe to call from multiple threads with *different* data | All solvers, algorithms |
| **Internally Synchronized** | Fully thread-safe (has internal locking) | `ThreadPool` |
| **Not Thread-Safe** | Requires external synchronization | File I/O classes |

---

## Thread-Safety by Layer

### BASE Layer (`mml/base/`)

| Class | Category | Notes |
|-------|----------|-------|
| `Vector<T>` | Thread-Safe (const) | Concurrent reads OK; writes need sync |
| `VectorN<T,N>` | Thread-Safe (const) | Fixed-size, same guarantees as Vector |
| `Matrix<T>` | Thread-Safe (const) | Concurrent reads OK; `Resize()` not safe |
| `MatrixNM<T,N,M>` | Thread-Safe (const) | Fixed-size matrix |
| `Complex<T>` | Thread-Safe | Immutable operations |

### CORE Layer (`mml/core/`)

| Class | Category | Notes |
|-------|----------|-------|
| `LUSolver` | Reentrant | Stateless; safe with different matrices |
| `QRSolver` | Reentrant | Stateless |
| `CholeskySolver` | Reentrant | Stateless |
| `SVDSolver` | Reentrant | Stateless |
| `NumDiff` | Reentrant | Pure numerical differentiation |
| `NumIntegration` | Reentrant | Pure numerical integration |
| `CoordTransf*` | Thread-Safe | Stateless coordinate transforms |
| `Derivation` | Reentrant | Pure function operations |
| `Integration` | Reentrant | Pure function operations |

### ALGORITHMS Layer (`mml/algorithms/`)

| Class | Category | Notes |
|-------|----------|-------|
| `RootFinder` | Reentrant | All root-finding methods |
| `BulirschStoerSolver` | Reentrant | Create one instance per thread |
| `RungeKuttaSolver` | Reentrant | Create one instance per thread |
| `DormandPrinceSolver` | Reentrant | Create one instance per thread |
| `EigenSolver` | Reentrant | Eigenvalue computations |
| `BVPSolver` | Reentrant | Boundary value problems |
| `InterpolatedFunction` | Thread-Safe (const) | Safe reads after construction |

### TOOLS Layer (`mml/tools/`)

| Class | Category | Notes |
|-------|----------|-------|
| `ThreadPool` | Internally Synchronized | Fully thread-safe |
| `Timer` | Internally Synchronized | Uses atomic counters |
| `Serializer` | Not Thread-Safe | File I/O; use one per thread |
| `DataLoader` | Not Thread-Safe | File I/O; use one per thread |
| `Visualizer` | Not Thread-Safe | Process spawning |
| `CSVExporter` | Not Thread-Safe | File I/O |

### Mathematical Functions

All pure mathematical functions in `MMLBase.h` are **thread-safe**:
- Trigonometric: `sin`, `cos`, `tan`, etc.
- Exponential: `exp`, `log`, `pow`
- Special functions: `gamma`, `beta`, `erf`
- Constants: `PI`, `E`, etc.

---

## Best Practices

### ✅ Safe Parallel Patterns

```cpp
// Pattern 1: Parallel ODE solving (one solver per thread)
#pragma omp parallel for
for (int i = 0; i < num_problems; i++) {
    BulirschStoerSolver solver;  // Thread-local instance
    auto result = solver.solve(problems[i]);
    results[i] = result;
}

// Pattern 2: Parallel matrix operations (read-only)
const Matrix<double> A = ...;
#pragma omp parallel for
for (int i = 0; i < n; i++) {
    double sum = 0;
    for (int j = 0; j < m; j++) {
        sum += A(i, j);  // Const access is thread-safe
    }
    row_sums[i] = sum;
}

// Pattern 3: Using ThreadPool
ThreadPool pool(4);
std::vector<std::future<double>> futures;
for (const auto& task : tasks) {
    futures.push_back(pool.enqueue([&task]() {
        return compute(task);
    }));
}
```

### ❌ Unsafe Patterns (Avoid)

```cpp
// WRONG: Shared mutable state
Matrix<double> shared_matrix(100, 100);
#pragma omp parallel for
for (int i = 0; i < 100; i++) {
    shared_matrix(i, i) = compute(i);  // Race condition!
}

// WRONG: Shared solver instance
BulirschStoerSolver shared_solver;  // Don't share!
#pragma omp parallel for
for (int i = 0; i < n; i++) {
    results[i] = shared_solver.solve(problems[i]);  // Race!
}

// WRONG: Concurrent file I/O
Serializer shared_serializer("output.csv");
#pragma omp parallel for
for (int i = 0; i < n; i++) {
    shared_serializer.write(data[i]);  // Race!
}
```

---

## Common Questions

### Q: Can I use OpenMP with MML?
**Yes.** Create solver instances inside the parallel region (thread-local). Read-only access to matrices/vectors is safe.

### Q: Is `Matrix::operator()` thread-safe?
**For reads (const), yes.** For writes, no—use synchronization or separate matrices per thread.

### Q: How do I parallelize ODE solving?
Create one solver instance per thread. The solvers maintain internal step history, so sharing would cause races.

### Q: Can multiple threads write to the same Serializer?
**No.** Create one Serializer per thread, or use external locking.

---

## Implementation Notes

- MML is a header-only library with no global mutable state
- Most classes are designed for single-threaded use with reentrant semantics
- `ThreadPool` is the only class with internal synchronization
- For maximum parallelism, prefer creating thread-local instances over sharing

---

## See Also

- [ThreadPool documentation](tools/ThreadPool.md)
- [Performance tips](QUICK_START.md#performance)
