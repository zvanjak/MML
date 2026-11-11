# Surface integration accuracy

`SurfaceIntegration` supports parametric midpoint flux and adaptive triangular
or rectangular mesh flux. Existing runnable calls are in
[docs_demo_integration_surface.cpp](../../src/docs_demos/core/docs_demo_integration_surface.cpp).

For adaptive mesh flux, `eps` is one absolute error budget for the complete
requested surface, including every face of a body. The algorithm refines the
largest available cell error and checks the sum of retained estimates against
that single budget. Passing a local tolerance on each child or face is not
sufficient for success.

An estimate compares a centroid contribution with the sum over four children.
The returned value, estimated error, deepest refinement actually used and actual
field evaluation count describe the same retained domain. A depth or numerical
stagnation limit yields an unusable result; the plain scalar wrapper throws
instead of discarding that status.

The existing default absolute tolerance is exposed as
`Precision::SurfaceIntegrationTolerance`, and the depth default as
`Precision::SurfaceIntegrationMaxDepth`. Sample-based estimates are not certified
mathematical error bounds.

## Parametric midpoint flux and scale

The parametric rule integrates the oriented scalar triple product
`F dot (r_u cross r_v)` directly. Binary-scaled products retain the field and
quadrature weights before reconstructing exponents; no unit normal or
absolute small-Jacobian cutoff is used. Triangular mesh contributions likewise
use one half of the oriented edge cross product, without normalizing it.
Reversing the parameterization or triangle winding reverses flux.

Range-stress regressions choose powers from the active `Real` exponent range,
keeping geometry, field values and parameter limits finite and nonzero while
the unweighted cross product deliberately overflows or underflows. Affine
finite-difference checks use a float-specific rounding allowance where needed;
double fixtures and bounds remain unchanged. For limits at the largest finite
float, the expected flux accounts for that endpoint's rounded mantissa rather
than treating it as an exact power of two. Orientation, callback counts and
result-status checks are independent of these numerical comparisons.

Requested parameter endpoints must be finite, strictly increasing and inside
the surface's declared rectangular domain. Subdivision counts must be positive,
and their product must fit the evaluation counter. Infinite declared domain
endpoints are permitted when the requested rectangle is finite.

Numerical partials use three distinct representable nodes inside the requested
rectangle. Centered quadratic stencils are preferred; one-sided quadratic
stencils handle midpoints rounded to a boundary. Actual node separations are
used rather than an assumed step. If three distinct nodes are unavailable,
`NumericalMethodError` is thrown. Nonfinite surface points, tangents, field
values or flux arithmetic also throw; callback exceptions propagate.

The midpoint result is `Completed`, not `Success`: `converged` and
`error_estimated` are false, `error_estimate` is `NoEstimate`, and `iterations`
retains the legacy patch count. `function_evaluations` counts actual field calls, one per patch,
including zero-area patches. Surface geometry is sampled five times per patch.
The triangular mesh rule also samples the field on exactly degenerate triangles;
zero area no longer bypasses callback validation.

Run `MML_DocsApp --geometric-integration-contracts` for checked examples:
`r(u,v)=(1e-8*u,1e-8*v,0)` with `Fz=1e16` has unit flux, as does
`r(u,v)=(u,v,0)` on `[0,1e-7]^2` with `Fz=1e14`.

These changes do not certify finite-difference accuracy. Coordinate cancellation,
unresolvable geometry and severe cancellation in flux sums can still produce
finite but inaccurate results with `Completed` status: this fixed rule checks no
accuracy criterion. Detected nonfinite arithmetic or collapsed stencil nodes
throw instead. The existing rectangular-mesh geometry constructor and normal
constraints are unchanged.