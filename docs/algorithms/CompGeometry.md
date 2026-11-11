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

All four return exactly `-1`, `0` or `1`: the sign of the corresponding
determinant of the represented finite binary coordinates, over `Real`'s
entire finite exponent range. This is **not** a determinant-value API:
a sign remains meaningful even when its nonzero magnitude cannot be
represented. Nonfinite coordinates in any argument throw `GeometryError`
before determinant or expansion arithmetic.

Adaptive evaluation retains proven floating-point filters and, for 3D
near-degeneracy cases, error-free differences, FMA residuals and expansions.
Their use is restricted to safe product ranges, including the smallest
significand bits needed by expansion products. Outside those ranges, fixed
binary integer accumulators sum exact signed products of the **original**
coordinates. No uniform rescaling can erase a small coordinate merely
because a different coordinate is large. The shared degree-four backend
serves `InCircle2D`; degree-three and degree-five determinants provide the
3D fallbacks. `Orientation2D` also shares its implementation with polygon
simplicity testing.

Expansion lengths and binary accumulator sizes are bounded at compile time
for the selected `Real` and stored in `std::array`; fallbacks require no
arbitrary-precision dependency or heap allocation. These guarantees assume
the normal binary floating-point environment with round-to-nearest,
subnormal preservation and no unsafe fast-math transformations.

Exact predicates cannot recover coordinate differences already rounded away
before the call. Inside/outside interpretations require a nondegenerate
defining triangle or tetrahedron and the specified orientation. A zero
determinant for degenerate defining points does not establish a unique
circumcircle or circumsphere. Robust signs do not make triangle areas,
volumes, circumcenters or interpolation values exact.

### Runnable mixed-scale predicate example

This complete example is compiled in
[docs_demo_predicate_contracts.cpp](../../src/docs_demos/algorithms/docs_demo_predicate_contracts.cpp).
Run `MML_DocsApp --predicate-contracts` to check it and the planar example below.

```cpp
#include <mml/algorithms/CompGeometry/RobustPredicates.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

void Docs_Demo_PredicateContracts() {
	using namespace MML;
	using MML::CompGeometry::RobustPredicates;

	const Real large = std::ldexp(Real(1), std::numeric_limits<Real>::max_exponent / 2);
	const Real tiny = std::ldexp(Real(1), std::numeric_limits<Real>::min_exponent + std::numeric_limits<Real>::digits);
	const Point3Cartesian a(0, 0, tiny), b(large, 0, 0), c(0, large, 0), origin(0, 0, 0);
	if (RobustPredicates::Orientation3D(a, b, c, origin) != 1 ||
		RobustPredicates::Orientation3D(a, c, b, origin) != -1)
		throw std::runtime_error("Mixed-scale orientation must retain the tiny coordinate");

	const Point3Cartesian on(large, large, tiny);
	const Point3Cartesian inside(large, large, std::nextafter(tiny, Real(0)));
	const Point3Cartesian outside(large, large, std::nextafter(tiny, Real(1)));
	if (RobustPredicates::InSphere3D(origin, c, b, a, on) != 0 ||
		RobustPredicates::InSphere3D(origin, c, b, a, inside) != 1 ||
		RobustPredicates::InSphere3D(origin, c, b, a, outside) != -1)
		throw std::runtime_error("In-sphere signs must distinguish adjacent represented offsets");

	std::cout << "Exact predicate range, cancellation and mixed-scale contracts verified.\n";
}
```

## Planar triangulation

`Triangulation::ComputeDelaunay` sorts distinct represented coordinates
lexicographically, constructs a triangulation by inserting visible hull edges,
then legalizes interior edges with robust Lawson flips. There is no finite
supertriangle whose removal could erase a skinny valid face. Sorting costs
O(n log n); hull scans and edge legalization can require quadratic work, with
ordered edge-map operations adding logarithmic overhead. This is not an
O(n log n) implementation.

The result preserves the original `points` array and its indices. Exact
duplicates use their first input index in the mesh; nearby distinct points
are not clustered by a tolerance. Fewer than three distinct points and
wholly collinear inputs return no triangles. Otherwise all distinct points,
including collinear hull-edge points, participate in nondegenerate CCW faces.
Cocircular diagonals are not flipped on ties; callers must not depend on a
particular diagonal or triangle ordering.

For triangle `t`, `neighbors[t][e]` is the neighbor across the directed edge
`triangles[t][e]` to `triangles[t][(e+1)%3]`, or `-1` on the hull. Interior
neighbors are reciprocal and use the opposite direction of the shared edge.
No auxiliary vertex indices appear in the result.

`LocatePoint` uses exact orientation signs and includes triangle boundaries.
`Interpolate` uses signed subtriangle determinants with separate binary
exponents, not an absolute area cutoff. Determinant magnitudes and interpolated
values remain floating-point approximations, not exact real arithmetic.
Constant samples are reproduced exactly. The caller's default is returned
only for a point outside the triangulated region (including an empty mesh).

**Invalid-input behavior:** nonfinite coordinates or sample values throw
`GeometryError`. A sample-count mismatch now throws `ArgumentError`, rather
than returning the outside-point default. A nonzero interpolation result that
cannot be represented also throws `GeometryError`. The public result vectors
must not be edited into inconsistent indices or adjacency.

`EarClipTriangulationDetailed` accepts simple polygons of either winding.
Winding is classified by an exact local turn at a lexicographic extreme,
avoiding cancellation in a globally translated shoelace sum. Returned valid
faces are CCW, including a three-vertex polygon. It is not a complete
self-intersection validator: callers must supply a simple polygon. If no ear
can be clipped, the existing diagnostic fan fallback remains
`InvalidInput`, with `used_fan_fallback`, `unclipped_vertices` and an error
message. The convenience `EarClipTriangulation` throws for that fallback;
fallback triangles must not be treated as a valid polygon triangulation.

These guarantees apply to coordinates as represented by `Real`; translating
or scaling before the call can itself round distinct points together or make
them nonfinite. Triangle area and Voronoi circumcenter measurements have
separate numerical contracts and are not made exact by these topology fixes.

### Circle and Voronoi construction contracts

`Circles::MinimumEnclosingCircle` and `Voronoi::ComputeVoronoi` validate finite
input coordinates before handling empty or singleton inputs. Circumcenters use
local coordinate differences and exponent-scaled arithmetic; collinearity is
decided by the existing exact determinant of represented input coordinates,
not an absolute area threshold. Tiny noncollinear triangles therefore remain
noncollinear, and translation does not introduce global squared-coordinate
cancellation. Center and radius measurements are rounded, not exact predicates.

The minimum-circle algorithm retains Welzl's randomized expected-O(n) method,
empty input's zero circle at the origin, and singleton input's zero-radius
circle at that point. Diameter midpoints avoid coordinate-sum overflow, and
half-distances are formed before converting to `Real`. Boundary containment
uses a relative roundoff allowance, with no unit-scale absolute padding.
A final pass encloses the inputs around the rounded center; tiny roundoff
adjustments to the minimum radius are possible. A center or radius that cannot
be represented raises `GeometryError`, rather than returning NaN or infinity.

Voronoi results preserve original site indices. Exact duplicates share the
first input occurrence's cell; edge site indices use that representative.
Distinct nearby sites are not clustered by tolerance. For `k` distinct
collinear sites, the result has `k-1` perpendicular bisectors between
lexicographically adjacent sites, no finite Voronoi vertices, and all cells
unbounded. Zero or one distinct site has no edges. An edge with `v1 == v2 == -1`
denotes a full line through the midpoint of `sites[site1]` and `sites[site2]`;
its unit `direction` spans that line. Other infinite edges are rays with a
finite vertex and an outward unit direction. Directions use scaled local
differences and exact orientation signs.

A noncollinear mesh unexpectedly lacking faces, a degenerate Delaunay face,
or an unrepresentable circumcenter raises `GeometryError`: no centroid is
substituted and no incomplete empty diagram is reported as successful.
These construction guarantees do not change the separate numerical contracts
of `Triangle2D` measurements.

`FindNearestSite` compares exponent-scaled distances without materializing raw
squares or an overflowing distance in `Real`. Its linear search returns the
first input index on a rounded-distance tie, including exact duplicate sites,
and `-1` for an empty diagram. Nonfinite query coordinates raise `GeometryError`
even for an empty diagram. Distances are rounded measurements, not exact
distance-order predicates for arbitrarily close ties.

`GetClippedCell` intersects the requested rectangle with every distinct-site
nearest-site halfplane, for both bounded and unbounded cells. It does not expand
the rectangle or depend on rounded circumcenters. Bisectors use local differences
and exponent-scaled products, not global coordinate sums, squared `Real` normals
or an absolute separation cutoff. Only exact duplicates omit a halfplane.
Signed halfplane evaluations within `Precision::RoundoffZeroThreshold` times
the sum of their absolute local terms are treated as boundary zero. The same
classification controls vertex retention and edge crossings; there is no
unit-scale absolute padding. Intersections retain scaled weights until the final
coordinate conversion, so an underflowing interpolation parameter does not erase
a representable vertex. Rounded coordinates are constrained to their source
segment's coordinate ranges. These are rounded polygons, not exact constructions;
near a bisector, roundoff and coordinate representability limit accuracy.

Clipped polygons retain CCW order and omit consecutive exact duplicate vertices.
Zero-area bounds, empty intersections, and intersections represented only by a
point or line return an empty polygon. Valid bounds with an invalid site index
also return empty. Bounds must be finite and ordered (checked before index/empty
shortcuts); otherwise `GeometryError` is raised. Exact duplicate sites receive
the same clipped cell. These query helpers do not change diagram topology.

```cpp
const Real scale = Real{1e-6};
const auto circle = Circles::MinimumEnclosingCircle(std::vector<Point2Cartesian>{
    {0, 0}, {2 * scale, 0}, {scale, std::sqrt(Real{3}) * scale}});
const Real offset = Real{1e9};
const auto translated = Voronoi::ComputeVoronoi({
    {offset, offset}, {offset + 2, offset}, {offset, offset + 2}});
const auto collinear = Voronoi::ComputeVoronoi({{2, 0}, {0, 0}, {1, 0}});

const int exponent = 3 * std::numeric_limits<Real>::max_exponent / 4;
const Real queryScale = std::scalbn(Real{1}, exponent);
const auto farSites = Voronoi::ComputeVoronoi({{0, 0}, {queryScale, 0}});
const int nearest = farSites.FindNearestSite({2 * queryScale, 0}); // 1
const Real spacing = std::scalbn(Real{1}, -exponent);
const auto tinySites = Voronoi::ComputeVoronoi({{0, 0}, {spacing, 0}, {2 * spacing, 0}});
const auto middle = tinySites.GetClippedCell(1, -spacing, -spacing, 3 * spacing, spacing);
// Four vertices; normalized x coordinates are 0.5 and 1.5, y coordinates are -1 and 1.
```

Runnable counterpart: `Docs_Demo_CircleVoronoiContracts()` in
[docs_demo_computational_geometry.cpp](../../src/docs_demos/algorithms/docs_demo_computational_geometry.cpp);
run `MML_DocsApp --circle-voronoi-contracts`.

### Runnable topology, winding and interpolation example

The following complete example is compiled in
[docs_demo_triangulation_contracts.cpp](../../src/docs_demos/algorithms/docs_demo_triangulation_contracts.cpp).
Build `MML_DocsApp` and run it with `--triangulation-contracts`.

```cpp
#include <mml/algorithms/CompGeometry/Triangulation.h>
#include <iostream>
#include <stdexcept>

void Docs_Demo_TriangulationContracts() {
	using namespace MML;
	using namespace MML::CompGeometry;

	const auto skinny = Triangulation::ComputeDelaunay({{0, 0}, {1, 0}, {Real(0.5), Real(1e-6)}});
	if (skinny.NumTriangles() != 1)
		throw std::runtime_error("A noncollinear three-point set must have one face");

	const Real side = Real(1e-6);
	const auto tiny = Triangulation::ComputeDelaunay({{0, 0}, {side, 0}, {0, side}});
	if (tiny.Interpolate({side / 4, side / 4}, {1, 1, 1}, -77) != 1 ||
		tiny.Interpolate({2 * side, 2 * side}, {1, 1, 1}, -77) != -77)
		throw std::runtime_error("Interpolation must distinguish tiny interior points from outside points");

	const Real offset = Real(8192);
	const Polygon2D clockwise({{offset, offset}, {offset, offset + 1},
		{offset + 1, offset + 1}, {offset + 1, offset}});
	const auto clipped = Triangulation::EarClipTriangulationDetailed(clockwise);
	if (!clipped.IsSuccess() || clipped.used_fan_fallback || clipped.triangles.size() != 2)
		throw std::runtime_error("Translation must not change polygon winding or trigger fallback");

	std::cout << "Planar triangulation topology, winding and interpolation contracts verified.\n";
}
```

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