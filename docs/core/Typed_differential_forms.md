# Typed Differential Forms

Typed differential forms add semantic wrappers around existing `VectorN`, `MatrixNM`, and tensor machinery. They are intended for geometric code where variance, frame, metric, and orientation matter. The untyped `VectorN` APIs remain the general-purpose numerical layer.

## Quick Reference

| Concept | Type or function | Meaning |
|---|---|---|
| Coordinate frame | `Cartesian3`, `Polar2`, `Spherical3` | Compile-time tag for components |
| Tangent vector | `TangentVector<N, Frame>` | Contravariant rank-1 object, velocity/displacement-like |
| Covector | `Covector<N, Frame>` | Covariant rank-1 object, gradient/one-form-like |
| Point | `Point<Real, N, Frame>` | Coordinates of a location in a frame |
| Differential form | `DifferentialForm<Scalar, N, K, Frame>` | Full-storage covariant alternating `K`-form |
| One-form | `Form1<N, Frame>` / `OneForm<N, Frame>` | Degree-1 `DifferentialForm`; bridges with `Covector` |
| Wedge product | `Wedge(a, b)` / `wedge(a, b)` | Exterior product with degree checks |
| Metric | `Metric<N, Frame>` | Positive-definite metric value object |
| Musical maps | `Flat` / `flat`, `Sharp` / `sharp`, `Inner` / `inner` | Lower and raise rank-1 components; compute metric inner product |
| Hodge star | `hodge_star(form, metric, orientation)` | Maps `K`-forms to `(N-K)`-forms |
| Typed cross product | `cross(a, b, metric, orientation)` | 3D cross product derived from forms |
| Scalar fields | `IScalarField<N, Frame>` | Typed scalar field evaluated at typed points |
| Exterior derivative | `exterior_derivative(field)` | Scalar field derivative as a one-form field |
| Form derivative | `exterior_derivative(formField)` | K-form field derivative as a `(K+1)`-form field |
| Metric gradient | `gradient(field, metric)` | Vector gradient as `sharp(df)` |
| Directional derivative | `DirectionalDerivative(field, point, direction)` | Convenience for `df(direction)` |
| Function field adapters | `ScalarFunctionFieldAdapter`, `VectorFunctionFieldAdapter` | Bridge existing scalar/vector function objects into typed field interfaces |
| Typed Cartesian field ops | `GradientCart`, `DivergenceCart`, `CurlCart`, `LaplacianCart` | Typed wrappers over Cartesian field operations |
| Coordinate maps | `CoordTransfCoordinateMap<...>`, common map aliases | Typed adapter over existing `CoordTransf` Jacobians |

## Core Rule

The metric-free derivative of a scalar field is a one-form. A vector gradient is a metric-dependent object:

```cpp
OneForm<3, Cartesian3> df = exterior_derivative(field)(point);
Real along_v = df(direction);

Metric<3, Cartesian3> metric = Metric<3, Cartesian3>::Euclidean();
TangentVector<3, Cartesian3> grad = gradient(field, metric)(point);
```

This is why `DirectionalDerivative` only needs a field, point, and tangent vector, while `gradient` also takes a metric.

The same metric-free exterior derivative is available for typed form fields below top degree:

```cpp
IFormField<2, 1, Cartesian2>& alpha = /* P dx + Q dy */;
Form2<2, Cartesian2> dAlpha = exterior_derivative(alpha)(point);
```

For a one-form `alpha = P dx + Q dy`, the component `dAlpha.Component(0, 1)` is `dQ/dx - dP/dy`. Top-degree form fields do not expose an exterior derivative overload.

These APIs use C++20 constraints at the overload boundary: forms evaluate only on the correct number of same-frame tangent vectors, wedge products are available only when `P + Q <= N`, the typed cross product is available only in 3D, and exterior derivatives of form fields are available only below top degree.

## Rank-1 Objects

Use `TangentVector` for contravariant components and `Covector` for covariant components. They share storage behavior with `VectorN`, but the type system rejects mixing frames or variance accidentally.

```cpp
TangentVector<3, Cartesian3> velocity{1.0, 0.0, 0.0};
Covector<3, Cartesian3> df{2.0, 2.0, -1.0};

Real directional = df(velocity);
```

Convert between `Covector` and degree-1 forms explicitly:

```cpp
OneForm<3, Cartesian3> oneForm = ToForm(df);
Covector<3, Cartesian3> covector = ToCovector(oneForm);
```

## DifferentialForm Class

The primary form type is `DifferentialForm<Scalar, N, K, Frame>`, where `K` is the form degree. It stores all `N^K` covariant components and enforces `0 <= K <= N` at compile time.

Common aliases cover the usual real-valued forms:

```cpp
ScalarForm<N, Frame>          // DifferentialForm<Real, N, 0, Frame>
Form1<N, Frame>               // DifferentialForm<Real, N, 1, Frame>
OneForm<N, Frame>             // readability alias for Form1<N, Frame>
Form2<N, Frame>               // DifferentialForm<Real, N, 2, Frame>
Form3<N, Frame>               // DifferentialForm<Real, N, 3, Frame>
VolumeForm<N, Frame>          // DifferentialForm<Real, N, N, Frame>
```

Components use one index per degree. A two-form in 3D therefore has components such as `omega.Component(0, 1)` and `omega.Component(2, 0)`.

`Component(...)` and `ComponentAt(...)` are raw component accessors. They do not automatically fill sign-related permutations. Use `SetAlternatingComponent(...)`, `IsAlternating()`, and `AlternatingPart()` when you want explicit alternating-form semantics.

```cpp
OneForm<3, Cartesian3> dx = BasisOneForm<0, 3, Cartesian3>();
OneForm<3, Cartesian3> dy = BasisOneForm<1, 3, Cartesian3>();

Form2<3, Cartesian3> area = Wedge(dx, dy);
Real xyComponent = area.Component(0, 1);

TangentVector<3, Cartesian3> ex{1.0, 0.0, 0.0};
TangentVector<3, Cartesian3> ey{0.0, 1.0, 0.0};
Real orientedArea = area(ex, ey);
```

Use `DifferentialForm` directly when a generic algorithm needs arbitrary degree `K`; use the aliases when the degree is part of the API contract.

## Alternating Form Helpers

`DifferentialForm` keeps full `N^K` storage for simple indexing and generic algorithms. Alternating-form invariants are explicit helpers:

```cpp
Form2<3, Cartesian3> omega;
omega.SetAlternatingComponent(REAL(4.0), 0, 1);

bool ok = omega.IsAlternating();        // true
Real xy = omega.Component(0, 1);        //  4
Real yx = omega.Component(1, 0);        // -4

Form2<3, Cartesian3> repaired = raw.AlternatingPart();
```

Use `SetAlternatingComponent` when constructing forms by independent components. Use `AlternatingPart` when a raw covariant component table should be projected to its alternating part.

## Exterior Derivative of Form Fields

`exterior_derivative` works in two layers:

- `IScalarField<N, Frame>` -> `IFormField<N, 1, Frame>`
- `IFormField<N, K, Frame>` -> `IFormField<N, K + 1, Frame>` when `K < N`

The implementation differentiates each component numerically and applies the alternating exterior derivative formula. It does not require a metric.

For example, in 2D:

```cpp
class Alpha : public IFormField<2, 1, Cartesian2>
{
public:
	Form1<2, Cartesian2> operator()(const Point<Real, 2, Cartesian2>& p) const override
	{
		Form1<2, Cartesian2> alpha;
		alpha.Component(0) = p[0] * p[1];       // P(x,y)
		alpha.Component(1) = p[0] * p[0] + p[1]; // Q(x,y)
		return alpha;
	}
};

Alpha alpha;
Form2<2, Cartesian2> dAlpha = exterior_derivative(alpha)(point);
```

The identity `d(d f) = 0` is also supported numerically within derivative tolerance.

## Typed Cartesian Field Operations

For Cartesian fields, typed wrappers bridge to the existing `FieldOperations` layer while preserving `Point`, `TangentVector`, and frame tags at the API boundary:

```cpp
Point<Real, 3, Cartesian3> point{1.0, 2.0, 3.0};

TangentVector<3, Cartesian3> grad = GradientCart(scalarField, point);
Real laplacian = LaplacianCart(scalarField, point);

Real div = DivergenceCart(vectorField, point);
TangentVector<3, Cartesian3> curl = CurlCart(vectorField, point);
```

These wrappers are convenience APIs for Cartesian coordinate calculations. For metric-aware geometric gradients, prefer `exterior_derivative(field)` plus `gradient(field, metric)`.

## Metrics and Hodge Star

`Metric` stores covariant components and their inverse. `Flat` / `flat` lowers a tangent vector to a covector; `Sharp` / `sharp` raises a covector to a tangent vector. `Inner` / `inner` computes the metric inner product.

```cpp
Metric<3, Cartesian3> metric = Metric<3, Cartesian3>::Euclidean();
TangentVector<3, Cartesian3> v{2.0, -1.0, 3.0};

OneForm<3, Cartesian3> velocityFlat = ToForm(flat(metric, v));
Form2<3, Cartesian3> flux = hodge_star(velocityFlat, metric, Orientation::Positive);
```

The Hodge star requires both metric and orientation. That is deliberate: changing either changes the meaning of flux and normals.

## Coordinate Maps

Typed coordinate maps adapt existing `CoordTransf` implementations without replacing them. The adapter uses the existing convention `J(row=target, col=source)`.

```cpp
CoordTransfPolarToCartesian2D transform;
PolarToCartesian2DMap map = MakePolarToCartesian2DMap(transform);

Point<Real, 2, Polar2> polar{2.0, Constants::PI / 2.0};
Point<Real, 2, Cartesian2> cart = map_point(map, polar);

TangentVector<2, Polar2> angular{0.0, 1.0};
TangentVector<2, Cartesian2> pushed = push_forward(map, angular, polar);

Covector<2, Cartesian2> dx{1.0, 0.0};
Covector<2, Polar2> pulled = pull_back(map, dx, polar);
```

Tangent vectors are pushed forward by `J v`; covectors and forms are pulled back by the dual action.

Available convenience aliases include:

```cpp
PolarToCartesian2DMap
CartesianToPolar2DMap
CylindricalToCartesian3DMap
CartesianToCylindrical3DMap
SphericalToCartesian3DMap
CartesianToSpherical3DMap
IdentityCoordinateMap<Frame, N>
```

The Jacobian convention is always `J(row=target coordinate, column=source coordinate)`, so for `x' = map(x)`, `J(i,j) = partial x'_i / partial x_j`.

## Migration From VectorN

Use typed forms where the geometric category is part of the contract:

| Existing style | Typed style |
|---|---|
| `VectorN<Real, N>` point coordinates | `Point<Real, N, Frame>` |
| `VectorN<Real, N>` velocity/displacement | `TangentVector<N, Frame>` |
| `VectorN<Real, N>` gradient components | `Covector<N, Frame>` or `Form1<N, Frame>` |
| `FieldOps::Gradient(...)` | `exterior_derivative(field)` plus `gradient(field, metric)` |
| Manual cross product | `cross(a, b, metric, orientation)` |
| Manual coordinate Jacobian multiplication | `push_forward(map, v, point)` / `pull_back(map, form, point)` |

Keep `VectorN` when writing coordinate-free numerical containers, generic solvers, or performance-oriented kernels where no geometric promise is being made.

## Tensor and Relativity Boundary

Rank-2 and higher runtime-variance tensor machinery remains the right layer for curvature, Ricci/Einstein tensors, electromagnetic tensors, and general relativity workflows. Typed rank-1 forms complement that machinery by making scalar derivatives, flux forms, normals, and frame transformations explicit at API boundaries.

See also:

- [Vector_Forms_Improvement.md](../Vector_Forms_Improvement.md)
- [Tensor_geometry_and_relativity.md](../algorithms/Tensor_geometry_and_relativity.md)
- [Coordinate_transformations.md](Coordinate_transformations.md)
- [Metric_tensor.md](Metric_tensor.md)

## Runnable Examples

- `Example09_TypedFormsGradient`: shows `df` as a one-form, `df(v)`, `DirectionalDerivative`, and `grad=sharp(df)`.
- `Example10_TypedFormsFlux`: shows `velocity_flat`, Hodge-derived flux form, typed cross product, and surface area form.