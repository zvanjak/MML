# Computational Geometry

## Robust predicates

`RobustPredicates` provides exact signs for topology-changing decisions:

- `Orientation2D(a, b, c)`: positive for counter-clockwise order.
- `InCircle2D(a, b, c, d)`: positive when `d` is inside a counter-clockwise
	circumcircle.
- `Orientation3D(a, b, c, d)`: positive when `d` is below the oriented plane
	`(a, b, c)`.
- `InSphere3D(a, b, c, d, e)`: positive when `e` is inside the circumsphere of
	a positively oriented tetrahedron `(a, b, c, d)`.

Each predicate first evaluates the determinant with a proven floating-point
error bound. Near a degeneracy, it recomputes the determinant using
floating-point expansions and returns its exact sign for the supplied `Real`
coordinates. Expansion lengths are bounded at compile time and stored in
`std::array`; the fallback performs no heap allocation.

Exact predicates recover the exact sign of represented floating-point inputs.
They cannot recover coordinate differences that were already rounded away or
products that underflow before entering the expansion arithmetic.

## Convex hull in 3D

`ConvexHull3DComputer::Compute` uses an incremental hull algorithm. Exact
`Orientation3D` signs control coplanarity, face visibility, face winding, and
outside-point assignment. Floating-point distances are used only to rank
already classified candidates. The centroid of the initial tetrahedron remains
an interior witness when new horizon faces are oriented.

The resulting faces use outward winding. `Volume`, `SurfaceArea`, and
`Contains` remain floating-point measurements rather than exact predicates;
large translations or extremely thin hulls can therefore lose measurement
accuracy even when the hull topology is correct.