# Parallel Integration

Finite-interval 1D quadrature and uniform N-dimensional Monte Carlo now use
MML's own [base ThreadPool](../../mml/base/ThreadPool.h). No OpenMP, std::async,
external-pool ownership or work-stealing scheduler is introduced.

Include `mml/core/Integration.h`, or the focused
`mml/core/Integration/IntegrationParallel.h` and
`mml/core/Integration/MonteCarloParallel.h` headers.

## Partitioned 1D quadrature

```cpp
ParallelIntegrationConfig quadrature;
quadrature.thread_count = 4;
quadrature.partitions = 32;
quadrature.method = GAUSS10KRONROD21;
auto integral = IntegrateParallelDetailed(
    [](Real x) { return x * x; }, 0, 1, quadrature);
```

This splits the finite interval into indexed contiguous partitions, executes the
existing `IntegrateCore` dispatcher independently on each, then reduces values
in ascending index order using compensated summation. Reversed bounds work.
Zero-width intervals return exact zero with `Completed` and no evaluations.
Unrepresentable/collapsed partitions are rejected instead of being silently lost.

`ParallelIntegrationConfig` inherits `IntegrationConfig` and execution settings:

| Field | Meaning |
|---|---|
| `method` | TRAP, SIMPSON, ROMBERG, GAUSS10 or GAUSS10KRONROD21 |
| `tolerance` | Positive, partition-local quadrature tolerance |
| `thread_count` | Workers; zero selects hardware concurrency, or one if unavailable |
| `partitions` | Zero selects four partitions per requested worker |
| `min_parallel_tasks` | Batches below this threshold run inline |

Worker count is capped by task count. One worker or one partition uses the same
partition plan inline without constructing a pool. For a cheap integrand, explicitly
select one worker or raise `min_parallel_tasks`: interval length alone does not
predict callback cost. There is no universal automatic cost estimator.

`ParallelIntegrationResult` extends `IntegrationResult` with `partitions_used` and
`threads_used`. Function evaluations are counted exactly. `iterations` is the sum
of partition refinement counts, **not global refinement depth**.

Converted callback failures preserve elapsed time and the actual number of
executed calls after scheduled work has drained. The same rule applies to
parallel Monte Carlo: failed blocks do not report the configured sample count
as if it had all executed. A partial payload is unavailable until the algorithm
has a complete-domain estimate; failure does not manufacture a computed zero.

- All locally checked criteria met: `Success`, `converged=true`.
- Composite fixed Gauss10/GK21: `Completed`, `converged=false`.
- Embedded GK error estimates remain available despite the lack of a checked
  accuracy target. Gauss10 has `NoEstimate`, not a fictional zero error.
- Error estimates are summed in absolute value only if every partition supplies
  one. A failed partition retains a failure status and an indexed error message.

Gauss10/GK21 become **composite rules**: 32 GK21 partitions evaluate 672 points,
not one 21-point rule. Partition-local tolerances do not establish a global
relative bound, particularly under cancellation. Prioritized adaptive subdivision
with a coordinated global error budget is deferred.

`IntegrateParallel` is the plain counterpart: value on a usable result,
`ConvergenceError` otherwise. Both callables and `IRealFunction` work.

## Parallel Monte Carlo

```cpp
ParallelMonteCarloConfig sampling;
sampling.seed = 123;
sampling.num_samples = 20000;
sampling.thread_count = 4;
sampling.use_antithetic = true;
VectorN<Real, 2> lower{0, 0}, upper{1, 1};
auto estimate = IntegrateMonteCarloParallelDetailed(
    [](const VectorN<Real, 2>& x) { return x[0] + x[1]; },
    lower, upper, sampling);
```

This supports uniform box sampling and optional antithetic pairs. It does not
parallelize stratified, importance or quasi-random methods.

`ParallelMonteCarloConfig` inherits the serial statistical settings and common
execution settings. Its additional fields are:

| Field | Meaning |
|---|---|
| `block_size` | Independent draws per indexed block; pairs when antithetic |
| `min_parallel_samples` | Smaller evaluation budgets execute the same blocks inline |
| `exception_policy` | Propagate or convert exceptions into failure statuses |

The total evaluation budget is exact. At least two independent draws are required:
two ordinary evaluations, or an even budget of at least four antithetic evaluations.
One antithetic pair is one statistical sample but two function evaluations.
`samples_used` and `function_evaluations` count evaluations, not pairs.
Positive counts must fit the existing integer evaluation diagnostic.

Each indexed block owns an `mt19937`, initialized through `seed_seq` from the master
seed and block index. Blocks do not consume the calling thread's `Random` engine.
Results expose `seed_used`, `blocks_used` and `threads_used`. Seed zero requests
entropy. Explicitly use a nonzero seed for reproducible calls.

**For a fixed seed, block size, inputs and build/platform**, changing worker count
or taking the inline fallback preserves the random block streams and ordered
statistics reduction. Changing block size changes the sampling plan. This does not
promise the old serial RNG sequence or cross-standard-library bitwise identity.
Integrands that introduce their own randomness/state need their own reproducibility
contract.

Block statistics are merged using the Chan-Golub-LeVeque/Welford formula, including
between-block mean differences. Block variances and standard errors are **not**
averaged. Antithetic variance is computed from pair averages.

The shared serial/parallel stopping test is:

standard error <= max(absolute_error, target_error * abs(integral estimate)).

`target_error` is finite and positive; `absolute_error` is finite and nonnegative,
defaulting to zero. An estimate near zero no longer bypasses the uncertainty
test. `Success` agrees with `converged`; an exhausted sample budget that misses the
criterion returns `MaxIterationsExceeded` with its estimate and explanatory message.
The standard error is statistical, not a rigorous deterministic error bound.
Zero-volume domains return exact zero/`Completed` without sampling.

`IntegrateMonteCarloParallel` returns the value or throws on an unusable result.

## Ownership, failures and thread safety

Every call owns a private fixed-size FIFO pool. Tasks write disjoint indexed slots,
the batch drains completely, and task exceptions are rethrown by lowest task index
rather than nondeterministic completion order. `ConvertToStatus` reports invalid
input, numerical failures, resource limits or callback failures explicitly.
An enqueue failure also drains already-submitted work before referenced state dies.

The caller's `AlgorithmContext` is snapshotted into each task and restored afterward,
including inline execution. Thus refinement budgets do not accidentally revert to
worker-thread defaults.

**The shared integrand must be reentrant and thread-safe.** A const interface does
not prove this. Captured mutable state, caches, counters and nested mutable
integration helpers need synchronization or separate ownership. No parallel
2D/3D nested, improper, path or surface integration is added here.

ThreadPool itself moved from tools to base without a forwarding shim. Failed
worker startup stops and joins workers already created before rethrowing.
External reusable pools/task groups remain deferred: existing pool-wide waits and
exception queues are not batch-isolated, and same-pool nested waits can deadlock.

## Runnable evidence and benchmarking

Both examples above have checked counterparts in
[docs_demo_parallel_integration.cpp](../../src/docs_demos/core/docs_demo_parallel_integration.cpp).
Build `MML_DocsApp` and run `--parallel-integration` (no GUI).

[parallel_integration_benchmark.cpp](../../testing/parallel_integration_benchmark.cpp)
compares identical indexed workloads inline and with up to four workers, including
private-pool startup/destruction. It uses 400-term transcendental callbacks, 1,024
Gauss10 partitions and 65,536 Monte Carlo evaluations. It verifies exact equality of
the compared results and reports best-of-three wall-clock measurements.
Compile as a C++20 optimized standalone executable with the repository include root.
Speed ratios are machine/load-dependent evidence, never unit-test assertions.

Observed on the local Windows Release build (2026-10-16, four workers, best of
three, pool lifecycle included):

| Workload | Inline | Pool | Observed ratio |
|---|---:|---:|---:|
| Composite Gauss10 | 13.522 ms | 3.794 ms | 3.56x |
| Monte Carlo | 118.225 ms | 41.374 ms | 2.86x |

Both comparisons used identical task plans and verified bitwise-equal results.
These measurements describe intentionally expensive callbacks, not expected
speedups for cheap polynomials or every integration workload.
