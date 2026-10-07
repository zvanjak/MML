# Approximation and Curve Fitting

**Files**: `mml/base/ChebyshevApproximation.h`, `mml/algorithms/CurveFitting.h`

MML supports near-minimax Chebyshev approximation and SVD-based linear
least-squares fitting. Both APIs retain their direct convenience functions and
also expose diagnostics for selecting approximation complexity and diagnosing
ill-conditioned fits.

For tabulated-data method selection and the shared interpolation result/config
contract, see [Interpolated Functions](../base/Interpolated_functions.md).
Compilable examples live in
[`docs_demo_chebyshev.cpp`](../../src/docs_demos/algorithms/docs_demo_chebyshev.cpp),
[`docs_demo_curve_fitting.cpp`](../../src/docs_demos/algorithms/docs_demo_curve_fitting.cpp),
and [`docs_demo_interpolated_functions.cpp`](../../src/docs_demos/base/docs_demo_interpolated_functions.cpp).

## Adaptive Chebyshev Approximation

`ApproximateChebyshevDetailed` builds a fixed-size approximation and reports
the coefficients, sampled maximum and RMS errors, and the absolute coefficient
tail. `ApproximateChebyshevAdaptive` doubles the term count until both the
sampled maximum error and coefficient-tail estimate meet the requested
tolerance.

```cpp
#include <mml/base/ChebyshevApproximation.h>

ChebyshevApproximationConfig config;
config.tolerance = 1e-10;
config.initial_terms = 8;
config.max_terms = 128;
config.validation_samples = 257;

auto result = ApproximateChebyshevAdaptive(
    [](Real x) { return std::exp(x); }, -1.0, 1.0, config);

if (result.IsSuccess()) {
    Real value = result.approximation(0.25);
}
```

The result reports `terms_used`, `degree`, `samples_used`,
`coefficient_tail_error`, `max_error_estimate`, and `rms_error_estimate`. If
`max_terms` is reached first, the status is
`AlgorithmStatus::ToleranceUnachievable` and the best computed approximation
is still returned.

For plotting and sampled error analysis, use `Evaluate(points)`,
`UniformGrid(count)`, or `SampleUniform(count, points, values)`. Uniform grids
include both domain endpoints and require at least two points.

## Weighted Least Squares

`WeightedGeneralLinearLeastSquares` minimizes

$$
\sum_i w_i\left(y_i-\sum_j c_j\phi_j(x_i)\right)^2
$$

by solving the row-scaled system $\sqrt{W}Ac=\sqrt{W}y$ with SVD.
Weights must be finite and strictly positive. The result includes the weighted
residual norm, weighted MSE, weighted $R^2$, effective rank, condition number,
and `weight_sum`. `WeightedPolynomialFit` supplies the power basis, while
`WeightedGeneralLinearFitDetailed` also reports structured status and timing.

## Regularized Least Squares

`RidgeGeneralLinearLeastSquares` solves

$$
\min_c \lVert Ac-y\rVert_2^2+\lambda\lVert c\rVert_2^2.
$$

The constant basis coefficient is unregularized by default. Set
`regularize_constant` to `true` to penalize it as well. Lambda must be finite
and non-negative.

`TikhonovGeneralLinearLeastSquares` accepts a diagonal vector $d$ and solves
the augmented system

$$
\begin{bmatrix}A\\D\end{bmatrix}c=
\begin{bmatrix}y\\0\end{bmatrix},\qquad D=\operatorname{diag}(d).
$$

This permits coefficient-specific penalties. Results report
`regularization_parameter` and `coefficient_norm` in addition to the usual fit
statistics.

## Orthogonal Bases

The following factories produce degree-zero through `max_degree` functions in
the requested family:

- `MakeChebyshevFitBasis`
- `MakeLegendreFitBasis`
- `MakeHermiteFitBasis`
- `MakeLaguerreFitBasis`

They return the same `Vector<std::function<...>>` accepted by unweighted,
weighted, ridge, and Tikhonov fitting APIs.

```cpp
auto basis = MakeLegendreFitBasis(5);
auto fit = RidgeGeneralLinearLeastSquares(x, y, basis, 1e-6);
```

Chebyshev and Legendre bases use their natural interval $[-1,1]$. Hermite and
Laguerre functions are defined on their standard unbounded and non-negative
domains respectively; scale input data appropriately for numerical stability.