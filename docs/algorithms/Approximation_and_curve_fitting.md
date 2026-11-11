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

## Generalized Least Squares: Correlated Observation Errors

`GeneralizedLinearLeastSquares(x, y, covariance, basis)` accepts a **full
observation covariance matrix** C, not a precision matrix or a vector of weights.
It minimizes

$$
Q=(y-Ac)^T C^{-1}(y-Ac),\qquad A_{ij}=\phi_j(x_i).
$$

The implementation factors C = L L^T, whitens the design and observations by
forward triangular solves, then solves the transformed least-squares system
using SVD. It does not form C inverse or the normal equations. Covariance is
normalized internally before factorization so a global change of measurement
units does not trigger an absolute pivot threshold.

```cpp
#include <mml/algorithms/CurveFitting.h>
using namespace MML;

Vector<Real> x{0, 1, 2, 3};
Vector<Real> y{0, 2, 4, 3};
Matrix<Real> covariance(4, 4, {
    2, 1, 0, 0,
    1, 2, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1
});
Vector<std::function<Real(Real)>> basis{
    [](Real) { return Real(1); },
    [](Real value) { return value; }
};
auto fit = GeneralizedLinearLeastSquares(x, y, covariance, basis);
auto detailed = GeneralizedLinearFitDetailed(x, y, covariance, basis);
// In basis order: intercept = 7/8, slope = 1; Q = 21/8.
```

Run `MML_DocsApp --gls` for the checked counterpart, `Docs_Demo_GeneralizedLeastSquares`
in [docs_demo_curve_fitting.cpp](../../src/docs_demos/algorithms/docs_demo_curve_fitting.cpp).

**Input contract**

- There must be m nonempty, matching x/y observations and between 1 and m basis
  functions. Both callable and `Vector<const IRealFunction*>` bases are supported.
- C must be m-by-m, finite, symmetric and numerically positive definite.
  Symmetry uses `PrecisionValues<Scalar>::IsMatrixSymmetricTolerance` after
  normalizing by the largest covariance entry; accepted roundoff asymmetry is
  averaged. No diagonal jitter, pseudoinverse of C, or diagonal-only fallback is
  applied. Singular, indefinite and numerically unresolved covariances are rejected.
- Observations, callable handles/pointers and evaluated basis values are validated.
  These safety checks are mandatory even if `check_finite` is false.
- Rank-deficient **designs** are supported by the SVD minimum-norm solution.
  This does not mean singular **covariances** are supported.
- Inputs are not modified; each basis is evaluated once per observation.

**Result semantics**

The direct function returns `GeneralLinearFitResult<Scalar>`, with coefficients
in basis order and diagnostics for the whitened design. `effective_rank` is the
SVD rank; `condition_number` is infinite for a rank-deficient design.

- `residual_norm` = sqrt(Q), the Mahalanobis residual norm.
- `mean_squared_error` = Q / m. This is **not** Q divided by the residual degrees
  of freedom, and it is not the existing weighted fit's Q / sum(weights).
- Generalized mean = (1^T C^-1 y) / (1^T C^-1 1).
- `r_squared` = 1 - Q / ((y - mean * 1)^T C^-1 (y - mean * 1)).
  It may be negative for a model without an intercept. If centered variation is
  zero, R-squared and adjusted R-squared are undefined and reported as NaN.
- Adjusted R-squared uses effective rank r: 1 - (1 - R^2)(m - 1)/(m - r),
  and is NaN when m <= r.
- `weight_sum` stays zero: full covariance is not a diagonal weight vector.

Identity covariance reproduces ordinary least squares. Diagonal covariance with
C_ii = 1 / w_i reproduces weighted coefficients, residual norm and R-squared.
Multiplying covariance by a positive scalar leaves coefficients, rank, condition
number and R-squared unchanged; Q and MSE scale inversely with that scalar.

The Detailed API returns `CurveFittingResult<Scalar>` with timing, evaluation count
and `Completed` on a usable direct fit (not a claim that an accuracy target was
checked). With `estimate_error=false`, R-squared, adjusted R-squared and MSE are
NaN; coefficients, residual norm, rank and condition number are still available.
Default exception policy propagates errors. `ConvertToStatus` reports invalid
inputs/covariance as `InvalidInput`, unrepresentable numerical diagnostics as
`NumericalInstability`, and callback runtime failures as `AlgorithmSpecificFailure`
with their error messages. Nonconvergence of the underlying SVD propagates or is
converted by the existing curve-fitting exception policy.

This adds correlated-error linear fitting only. Existing diagonal-weighted,
ordinary, polynomial and regularized APIs retain their behavior. It does not
estimate covariance from residuals or implement nonlinear/feasible GLS.

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