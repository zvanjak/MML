# Stiff solver work-precision benchmarks

The benchmark suite measures Backward Euler, BDF2, and ROS2 against independently
verified endpoint values. Its purpose is to expose accuracy, work, and failure
behavior before and after implicit-solver changes, not to certify every configuration.

## Problems and references

| Problem | Initial state | Full endpoint | Source |
|---|---|---|---|
| Robertson | (1, 0, 0) | 10^11 | [Bari/CWI IVP test set](https://archimede.uniba.it/~testset/src/problems/rober.f) |
| Van der Pol, mu=1000 | (2, 0) | 2000 | [Bari/CWI IVP test set](https://archimede.uniba.it/~testset/src/problems/vdpol.f) |
| HIRES | (1, 0, 0, 0, 0, 0, 0, 0.0057) | 321.8122 | [Bari/CWI IVP test set](https://archimede.uniba.it/~testset/src/problems/hires.f) |

These classical stiff problems are also described in Hairer and Wanner,
*Solving Ordinary Differential Equations II*.

The published Van der Pol problem is in slow time with epsilon=10^-6. To compare
it with MML's unscaled equations at mu=1000, use physical time
t=1000*tau and physical velocity y[1]=v/1000. The published endpoint tau=2
therefore corresponds to t=2000, not t=2.

The short ("smoke") profile ends at t=1 for each problem. Its references were
computed independently using SciPy Radau and cross-checked using LSODA.
[generate_stiff_references.py](../../testing/generate_stiff_references.py)
reproduces that check and also checks the full endpoints against published values.
It requires NumPy and SciPy; the C++ tests and benchmark runner do not.

## Running and interpreting measurements

Build the Release target `MML_StiffBenchmarks`, then invoke the built executable
with a new output CSV filename, for example `results\stiff_work_precision_current.csv`.
Do not overwrite the committed historical baseline.
The runner is excluded from the default build. The ordinary `MML_Tests` target
includes the short reference, refinement, and work-counter tests under the
`[stiff-benchmark]` tag.

The runner records both profiles for all three methods at three settings (54 rows).
For fixed-step methods the initial step is the profile's duration divided by
100, 200, or 400. Full ROS2 runs start at duration/10000, with three progressively
tighter tolerances. This is intentionally a work-precision sweep, not an assertion
that fixed-step methods can accurately span Robertson's entire time interval with
a few hundred steps.

Each row records the requested endpoint, actual reached endpoint when available,
configuration, solver status, diagnostic message, maximum absolute endpoint error,
maximum component error scaled by tolerance*(1+abs(reference)), RHS and Jacobian
calls, accepted/rejected steps, and wall-clock milliseconds.
Work counters are measured through a wrapper around the problem; they do not
trust the legacy detailed APIs' estimated function counts. Timing is informational,
not a regression gate. Compare work and accuracy together, on the same profile.

Error estimates are NaN if a solver fails or does not reach the reference endpoint.
Failures are explicit CSV rows and stderr messages. The runner returns 1 when any
solver configuration fails, 2 on output/usage errors, and 0 only when all runs
succeed. A successful status alone does not imply adequate endpoint accuracy:
inspect the error columns. If a legacy solver discards its partial solution on
failure, the reached time is unavailable rather than inferred from its message.
Unavailable reached times are NaN; unavailable accepted/rejected step counts are
-1. Measured RHS/Jacobian calls remain available even in those failed runs.

BE/BDF2 enforce their configured work budget across accepted and rejected attempts
and retain partial trajectories on terminal status failures. ROS2 is configured
with a 50000-step capacity because the current
default 100000-step result allocation exceeds the matrix storage limit after
including the initial point. The current ROS2 driver also uses an internal
50000-step bound and a maximum step of 100 times its initial step; the benchmark
does not claim its other configuration bounds are enforced.

## Initial baseline

[stiff_work_precision_baseline.csv](../../results/stiff_work_precision_baseline.csv)
captures the pre-hardening results. Newton failures, singular solves, or exhausted
step budgets in full runs are evidence for the hardening epic, not silently
discarded measurements. Routine tests use shorter intervals so a full-suite run
does not become a long-running performance experiment.

## Newton-retry checkpoint

[stiff_work_precision_newton_retries.csv](../../results/stiff_work_precision_newton_retries.csv)
records the same 54 configurations after bounded BE/BDF2 Newton retries. There are
41 successful runs and 13 explicit failures, compared with 39 and 15 in the baseline.
Both formerly failing HIRES runs at duration/200 now reach the published endpoint.
Their maximum absolute errors are approximately 1.32e-5 (BE) and 2.89e-4 (BDF2).

The six full-period Van der Pol BE/BDF2 runs now recover from Newton failures and
advance farther, but exhaust the 50000-attempt budget at approximately t=815..867.
The retained smaller step is not grown again. Six full Robertson runs still
terminate on singular iteration matrices, and one Robertson ROS2 run still
exhausts its internal step budget. The benchmark therefore still returns 1.
These are explicit remaining limitations, not evidence that retries or all
full-profile accuracy targets have been solved.

## Factorization-reuse checkpoint

[stiff_work_precision_factor_reuse.csv](../../results/stiff_work_precision_factor_reuse.csv)
records step-local simplified Newton with LU refresh on slow convergence.
[stiff_work_precision_full_newton.csv](../../results/stiff_work_precision_full_newton.csv)
is the matched control: run the same executable with `--full-newton` after the
output filename to refresh the Jacobian and LU on every Newton correction.
The switch does not affect ROS2, which already shares one LU between its stages.
Both reports include the requested `reuse_jacobian` setting, successful
`lu_factorizations`, and completed `linear_solves`, including rejected attempts.
The historical reports retain their original schemas and measurements.

Representative matched work counts are:

| Problem/profile | Method/level | Fresh-Newton LU | Reused LU | Reduction |
|---|---|---:|---:|---:|
| Van der Pol smoke | BE / 0 | 104 | 100 | 3.85% |
| Van der Pol smoke | BDF2 / 0 | 103 | 100 | 2.91% |
| HIRES full | BE / 1 | 904 | 419 | 53.65% |
| HIRES full | BDF2 / 1 | 606 | 416 | 31.35% |

Van der Pol spends much of these configurations' work in nearly linear solves,
so its additional LU reduction is modest (0.28..1.11% on the failed full runs).
Avoiding Jacobians at already-converged residuals also reduces Jacobian counts;
that saving is distinct from reusing a factorization and must not be reported as
an equivalent LU saving. Simplified Newton may need more RHS evaluations and
back-substitutions than full Newton. Neither fewer LUs nor these timings promise
a wall-clock speedup for every problem size or RHS cost.

Across the 41 successful matched rows, the largest difference in the reported
maximum endpoint error is approximately 1.37e-11. Routine regressions also compare
the componentwise endpoints, constraints, and nonlinear residuals, including both
full-interval HIRES recovery runs with at least 25% fewer factorizations.

Both modes still record 13 failures and return 1. The six full Robertson runs
now get past the Gauss-Jordan singular-solve failures using scaled-pivot LU and
bounded retries, but exhaust their 50000-attempt budget at approximately t=727.
The six full Van der Pol runs and one full Robertson ROS2 run also remain
budget-limited. No new full-profile success or accuracy certification is claimed.
