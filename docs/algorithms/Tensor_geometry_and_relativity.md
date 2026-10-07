# Tensor Geometry and Relativity Algorithms

This note records the numerical conventions behind MML's tensor, differential-geometry, and relativity infrastructure, plus the local MPL electromagnetism extension built on it. It complements the API references in [../base/Tensors.md](../base/Tensors.md) and [../core/Metric_tensor.md](../core/Metric_tensor.md), and the roadmap documents [../Tensors_improvement.md](../Tensors_improvement.md) and [../Vector_Forms_Improvement.md](../Vector_Forms_Improvement.md).

## Scope

The current implementation covers two connected layers:

- Fixed-rank tensors `Tensor2<N>` through `Tensor5<N>` with runtime index variance.
- Metric-driven geometry in `MetricTensorField<N>`, including Christoffel symbols, Riemann/Ricci/Einstein curvature, geodesic integration, and parallel transport.

The companion vector/forms proposal covers the typed rank-1 layer: tangent vectors, covectors, one-forms, frames, wedge products, Hodge star, and explicit `flat`/`sharp` metric maps. The intended architecture is additive: forms make the user-facing rank-1 calculus safer, while the rank-2+ tensor engine remains a compact runtime-variance engine for curvature and physics.

## Conventions

### Index Order

Spacetime APIs use component order:

```text
(ct, x, y, z)
```

Schwarzschild examples use:

```text
(t, r, theta, phi)
```

The unit 2-sphere uses:

```text
(theta, phi)
```

### Metric Signature

Lorentzian spacetime helpers use the mostly-plus convention:

```text
eta_munu = diag(-1, +1, +1, +1)
```

Timelike four-velocities therefore satisfy `g(u,u) = -1` in natural units. The GR examples use this as a runtime sanity check for Schwarzschild geodesics.

### Riemann Tensor Convention

`MetricTensorField<N>` computes the mixed-index convention:

```text
R^rho_sigma_mu_nu = partial_mu Gamma^rho_sigma_nu
                   - partial_nu Gamma^rho_sigma_mu
                   + Gamma^rho_lambda_mu Gamma^lambda_sigma_nu
                   - Gamma^rho_lambda_nu Gamma^lambda_sigma_mu
```

The Ricci tensor contracts the first and third Riemann slots:

```text
R_sigma_nu = R^rho_sigma_rho_nu
```

and scalar curvature is:

```text
R = g^sigma_nu R_sigma_nu
```

With this convention, the unit 2-sphere metric `diag(1, sin^2(theta))` has `R^theta_phi_theta_phi = sin^2(theta)` and scalar curvature `R = 2`.

## Tensor Algebra Notes

Tensor classes store dimension and rank at compile time, but variance at runtime. That is deliberate:

- Curvature workflows constantly raise, lower, and contract indices.
- Runtime variance lets one contraction routine handle `T^i_j`, `T_i^j`, and higher-rank mixed tensors without an overload explosion.
- Compile-time frame and variance safety is still valuable at rank 1, where the forms proposal can prevent accidental tangent/covector/frame mixing.

Important operations and checks:

- Addition and subtraction require identical variance patterns.
- Contraction requires opposite variance on the contracted slots.
- `Tensor3<N>` may store Christoffel-like component arrays, but Christoffel symbols are not tensors and must not be transformed as tensor fields.
- The Levi-Civita symbol and metric-weighted Levi-Civita tensor are distinct objects. Use the tensor form when the metric determinant matters.

## Curvature Algorithm

The curvature pipeline is:

1. Evaluate covariant metric components `g_ij(q)`.
2. Invert `g_ij` to obtain `g^ij(q)`.
3. Numerically differentiate metric components to compute Christoffel symbols.
4. Numerically differentiate Christoffel symbols and combine quadratic connection terms for Riemann curvature.
5. Contract Riemann to Ricci, contract Ricci with inverse metric to scalar curvature, then form the Einstein tensor.

Reference validation geometries:

- Flat Cartesian, cylindrical, and spherical coordinates: all curvature quantities vanish even when Christoffel symbols are nonzero.
- Unit 2-sphere: constant positive intrinsic curvature.
- Schwarzschild exterior metric: nonzero Riemann curvature with numerically vanishing Ricci tensor, Ricci scalar, and Einstein tensor.

## Geodesic Algorithm

`GeodesicEquationSystem<N>` turns the second-order geodesic equation into a first-order ODE system with state layout:

```text
(q^0, ..., q^(N-1), v^0, ..., v^(N-1))
```

The right-hand side is:

```text
dq^i/dlambda = v^i
dv^i/dlambda = -Gamma^i_jk(q) v^j v^k
```

Use `IntegrateGeodesicFixedStep` for a direct fixed-step RK integration. For Lorentzian geodesics, monitor `g(v,v)` during validation. A timelike proper-time parametrized geodesic should preserve `g(v,v) = -1` up to numerical error.

## Electromagnetism Algorithm

The MPL flat-spacetime EM helpers in `src/book/mpl/Electromagnetism/Electromagnetism.h` use covariant four-potential components `A_mu` and compute:

```text
F_munu = partial_mu A_nu - partial_nu A_mu
```

Field extraction follows the mostly-plus convention used by `MetricTensorMinkowski`. The two invariants are:

```text
F_munu F^munu
epsilon^munurhosigma F_munu F_rhosigma
```

Liénard-Wiechert potentials solve the retarded-time equation for a moving point charge and then build the four-potential. The boosted Coulomb example validates the uniform-motion case by comparing fields extracted from `F_munu` against the analytic boosted Coulomb field.

## Numerical Accuracy

Most tensor-geometry operations use finite differences. Accuracy depends on scale, smoothness, and coordinate conditioning.

Practical guidance:

- Choose evaluation points away from coordinate singularities and metric degeneracies.
- Use problem scales near order one when possible, or tune derivative steps explicitly.
- Validate tensor identities by tolerance, not exact equality.
- Prefer invariant checks, such as scalar curvature, Ricci flatness, `g(u,u)`, and EM invariants, over isolated component checks when diagnosing numerical drift.

## Singularity Handling

Coordinate singularities and physical singularities must be treated differently.

- Spherical and cylindrical coordinate axes can be coordinate singularities even in flat space.
- The Schwarzschild horizon `r = r_s` is singular in Schwarzschild coordinates.
- The Schwarzschild physical singularity at `r = 0` is not a removable coordinate artifact.
- Liénard-Wiechert potentials are singular at the source worldline and can be ill-conditioned near null-cone tangencies.

Examples and tests should therefore pick points safely outside these regions. If an application must approach them, use a coordinate chart or formulation designed for that region instead of expecting finite-difference Christoffels to remain well-conditioned.

## Link to Differential Forms

The vector/forms roadmap supplies the missing semantic layer for differential forms and rank-1 objects:

- `df` should be a one-form, not an ordinary vector.
- `gradient(f)` should mean `sharp(df)` and should require a metric.
- Cross product should be derived from wedge product, Hodge star, metric, and orientation.
- EM can eventually expose `F` as a 2-form while still using the existing `Tensor2<4>` storage and invariant machinery underneath.

This division keeps numerical tensor algorithms compact while making the public calculus surface harder to misuse.