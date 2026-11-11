# Geometry - Points and Pure Geometric Shapes

**Files**: `mml/base/Geometry.h`, `mml/base/GeometrySpherical.h`

Lightweight point types for 2D/3D coordinate systems, plus **coordinate-free pure geometric shapes**.

## Table of Contents
- [2D Points](#2d-points)
- [3D Points](#3d-points)
- [Pure Geometric Shapes](#pure-geometric-shapes)
- [Examples](#examples)

---

## 2D Points

### Point2Cartesian

2D Cartesian point (x, y).

```cpp
// Construction
Point2Cartesian p1;           // Default: (0, 0)
Point2Cartesian p2(3.0, 4.0); // (3, 4)

// Access
Real x = p2.X();
Real y = p2.Y();

// Modification
p2.X() = 5.0;
p2.Y() = 12.0;

// Distance
Real dist = p1.Dist(p2);  // Distance between points
```

**Operations:**
```cpp
Point2Cartesian a(1, 2), b(3, 4);

// Arithmetic
auto sum = a + b;     // (4, 6)
auto diff = a - b;    // (-2, -2)
auto scaled = a * 2;  // (2, 4)
auto scaled2 = 3 * a; // (3, 6)
auto divided = a / 2; // (0.5, 1)

// In-place
a += b;
a -= b;
a *= 2.0;
a /= 2.0;

// Comparison
bool exact_match = (a == b);
bool not_equal = (a != b);
bool approx_equal = a.IsEqual(b, 1e-6);  // With tolerance
```

### Point2Polar

2D polar point (r, φ).

```cpp
// Construction
Point2Polar p1;                    // Default: (0, 0)
Point2Polar p2(5.0, Constants::PI / 4);  // r=5, φ=π/4

// Convert from Cartesian
Point2Cartesian cart(3, 4);
Point2Polar polar(cart);  // r ≈ 5, φ ≈ 0.927 rad

// Access
Real r = p2.R();
Real phi = p2.Phi();

// Convert to Cartesian
Point2Cartesian back = p2.TransfToCart();

// Distance (using law of cosines)
Real dist = p1.Dist(p2);
```

---

## 3D Points

### Point3Cartesian

3D Cartesian point (x, y, z).

```cpp
// Construction
Point3Cartesian p1;                // Default: (0, 0, 0)
Point3Cartesian p2(1.0, 2.0, 3.0); // (1, 2, 3)

// Access
Real x = p2.X();
Real y = p2.Y();
Real z = p2.Z();

// Modification
p2.X() = 4.0;
p2.Y() = 5.0;
p2.Z() = 6.0;

// Distance
Real dist = p1.Dist(p2);  // Euclidean distance
```

**Operations:**
```cpp
Point3Cartesian a(1, 2, 3), b(4, 5, 6);

// Arithmetic
auto sum = a + b;     // (5, 7, 9)
auto diff = a - b;    // (-3, -3, -3)
auto scaled = a * 2;  // (2, 4, 6)
auto scaled2 = 3 * a; // (3, 6, 9)
auto divided = a / 2; // (0.5, 1, 1.5)

// In-place
a += b;
a -= b;
a *= 2.0;
a /= 2.0;

// Comparison
bool exact_match = (a == b);
bool not_equal = (a != b);
bool approx_equal = a.IsEqual(b, 1e-6);
```

### Point3Spherical

3D spherical point (r, θ, φ).

```cpp
// Construction
Point3Spherical p(5.0, Constants::PI/2, 0.0);  // r, θ (theta), φ (phi)

// Access
Real r = p.R();
Real theta = p.Theta();  // Polar angle (from z-axis)
Real phi = p.Phi();      // Azimuthal angle

// Convert from/to Cartesian
Point3Cartesian cart(0, 0, 5);
Point3Spherical sph(cart);

Point3Cartesian back = sph.TransfToCart();
```

### Point3Cylindrical

3D cylindrical point (ρ, φ, z).

```cpp
// Construction
Point3Cylindrical p(3.0, Constants::PI/4, 2.0);  // ρ (rho), φ (phi), z

// Access
Real rho = p.Rho();  // Radial distance from z-axis
Real phi = p.Phi();  // Azimuthal angle
Real z = p.Z();      // Height

// Convert from/to Cartesian
Point3Cartesian cart(3, 0, 2);
Point3Cylindrical cyl(cart);

Point3Cartesian back = cyl.TransfToCart();
```

---

## Pure Geometric Shapes

Coordinate-free shapes defined by intrinsic parameters only (no position/orientation).

### 2D Shapes

| Class | Parameters | Key Methods |
|-------|------------|-------------|
| `Triangle` | sides a, b, c | `Area()`, `Perimeter()`, angles, medians, `Inradius()`, `Circumradius()`, `IsRight()`, `IsEquilateral()` |
| `Circle` | radius | `Area()`, `Circumference()`, `ArcLength(θ)`, `SectorArea(θ)`, `SegmentArea(θ)` |
| `Ellipse` | semi-axes a, b | `Area()`, `Circumference()`, `Eccentricity()`, `FocalDistance()` |
| `CircularSector` | radius, angle | `Area()`, `ArcLength()`, `ChordLength()` |
| `CircularSegment` | radius, angle | `Area()`, `Sagitta()`, `FromSagitta()` factory |
| `Annulus` | inner/outer radius | `Area()`, inner/outer circumferences |
| `RegularPolygon` | n sides, length | `Area()`, `Apothem()`, `Circumradius()`, named factories |
| `Rectangle` | width, height | `Area()`, `Diagonal()`, `AspectRatio()`, `GoldenRectangle()` |
| `Parallelogram` | sides a, b, angle | `Area()`, diagonals, heights |
| `Rhombus` | diagonals d1, d2 | `Area()`, angles, `Inradius()` |
| `Trapezoid` | parallel sides, height | `Area()`, `Median()`, `Isosceles()` factory |

### 3D Shapes

| Class | Parameters | Key Methods |
|-------|------------|-------------|
| `SphereGeom` | radius | `Volume()`, `SurfaceArea()`, `CapVolume(h)`, `GreatCircle()` |
| `CylinderGeom` | radius, height | `Volume()`, `LateralArea()`, `SurfaceArea()` |
| `ConeGeom` | radius, height | `Volume()`, `SlantHeight()`, `ApexAngle()` |
| `Frustum` | two radii, height | `Volume()`, `SlantHeight()`, `SurfaceArea()` |
| `Tetrahedron` | 6 edge lengths | `Volume()` (scaled Gram determinant, equivalent to Cayley-Menger), `SurfaceArea()`, face triangles |
| `Spheroid` | equatorial, polar radii | `Volume()`, `Flattening()`, `Eccentricity()`, `Earth()` |
| `TorusGeom` | major, minor radius | `Volume()`, `SurfaceArea()`, `IsRingTorus()` |

### Quick Examples

```cpp
// Triangle from sides (Heron's formula)
Triangle t(3, 4, 5);
Real area = t.Area();           // 6.0
bool right = t.IsRight();       // true
Real inr = t.Inradius();        // 1.0

// Circle
Circle c(5.0);
Real circ = c.Circumference();  // 2πr ≈ 31.42
Real sector = c.SectorArea(Constants::PI/2);  // quarter circle

// Ellipse
Ellipse e(5, 3);
Real ecc = e.Eccentricity();    // 0.8

// Regular polygon
auto hex = RegularPolygon::Hexagon(10.0);
Real apothem = hex.Apothem();

// 3D shapes
SphereGeom s(1.0);
Real vol = s.Volume();          // 4π/3 ≈ 4.19

Tetrahedron tet = Tetrahedron::Regular(1.0);
Real tetVol = tet.Volume();     // Regular tetrahedron volume
```

> **Note**: 3D shape classes use `Geom` suffix (e.g., `SphereGeom`) to avoid conflicts with parametric surface classes in `Surfaces.h`.

### Shape and Mesh Contracts

- Function-bounded solids preserve their existing constructors, but
  `getBoundingBox()` and `getBoundingSphere()` require a conservative box
  supplied through `SetBoundingBox()`. The caller guarantees enclosure of the
  **entire** boundary functions; sampling cannot establish that guarantee.
  Without it, these queries throw `GeometryError`. The sphere encloses the
  supplied box, using its midpoint rather than an estimated centroid and
  outward-rounded radius arithmetic. Known changes to the public x interval
  or boundary callback pointers require a new box; callers must also maintain
  the guarantee when callback state changes externally. `ClearBoundingBox()`
  removes the supplied guarantee.
- `EstimateBoundingBox()` and `EstimateBoundingSphere()` are explicitly
  heuristic sampled extrema and can miss narrow peaks. `Volume()`,
  `SurfaceArea()` and `getCenter()` remain numerical estimates, forwarding to
  `EstimateVolume()`, `EstimateSurfaceArea()` and `EstimateCenter()`.
  The estimate methods accept `BoundaryGeometrySettings`: a positive even
  `integrationIntervals` count (default 64) and a finite `derivativeStep` in
  (0, 1/2], measured in normalized [0,1] parameters (default
  `Precision::CbrtEps`). These are resolution controls, **not** certified
  error bounds or convergence claims; fixed-grid quadrature can alias
  oscillatory callbacks even when all evaluations are finite.
- Function-bounded surface area differentiates the domain-preserving
  parameterization, including moving y boundaries, rather than sampling
  Cartesian partial derivatives outside the body domain. Three-point
  stencils are central in the interior and one-sided at endpoints; all
  physical samples stay within the closed callback domain. Callbacks must be
  deterministic and finite there, with ordered lower/upper boundaries and
  sufficient smoothness for meaningful numerical derivatives. Null callbacks,
  invalid settings, nonfinite/reversed sampled bounds, unresolved parameter
  stencils and nonfinite arithmetic throw `GeometryError`. A zero estimated
  volume has no centroid; small positive volumes are not treated as zero by an
  absolute geometry tolerance.
- `Polygon2D::IsSimple()` tests strict topology of the supplied finite
  coordinates, without an epsilon-expanded boundary. The last vertex closes
  back to the first implicitly: repeating the first vertex at the end is a
  zero-length edge and is not simple. Other zero-length edges, repeated
  vertices, crossings, non-adjacent touches, collinear overlaps and adjacent
  retracing are also rejected. Adjacent edges may share only their common
  endpoint; harmless collinear subdivisions and valid concavity remain simple.
  Shared exact-sign orientation arithmetic handles cancellation, product
  overflow/underflow and mixed exponent ranges without changing vertex storage.
  Nonfinite coordinates throw `GeometryError`. This predicate does not change
  containment, convexity or the existing convex-only fan triangulation contract.
- Cartesian point distances and Cartesian-to-polar/spherical/cylindrical
  constructors use scaled `hypot` norms. Curvilinear point distances combine
  native radial differences with scaled half-angle chord components rather
  than subtracting large squared terms. Signed radii and noncanonical finite
  angles remain accepted; spherical zero/pi poles and canonical angular wrap
  boundaries retain their geometric identities. Distance/radius operations
  reject nonfinite coordinates and unrepresentable required results with
  `GeometryError`, rather than returning a spurious coincidence or nonfinite
  distance. These changes do not redesign general vector norms or coordinate
  transformation APIs.
  Coordinates are interpreted after conversion to the configured `Real` type:
  a radial increment smaller than its representable spacing can round to the
  original radius, giving genuinely coincident stored points. No distance
  algorithm can recover that lost input difference. Extreme-coordinate inputs
  and near-antipodal offsets must themselves be finite and distinct in the active
  precision; this is separate from preserving a representable result when
  intermediate squares or angular products would overflow or underflow.
- Great-circle distances evaluate the near and complementary half-angle
  chords independently and use `atan2`, preserving both nearby and near-antipodal
  separations. Geographic latitudes must be finite and in [-90, 90]; finite
  longitudes are periodic, including outside [-180, 180]. The four-argument
  overload returns radians; the five-argument overload returns arc length in
  the supplied radius units. Radius must be finite and nonnegative; zero is
  allowed. The physical-distance overload scales tiny chord components before
  they underflow, even when the angular result itself is unrepresentable.
  Invalid inputs and unrepresentable required distances throw `GeometryError`.
  `SphericalFromLatLong` uses the same geographic validation and continues to
  store the supplied longitude in radians without canonicalizing it.
- `Triangle2D` computes side lengths from its current vertices rather than a
  mutable cache. Writes through retained vertex or coordinate references are
  visible to subsequent side-derived queries, including on copied/moved
  triangles. Read-only queries do not write shared cache state; concurrent
  vertex mutation still requires caller synchronization.
  `Area()` is the absolute coordinate-determinant `SignedArea()`, so a thin
  triangle does not lose its height merely because its rounded side lengths
  satisfy a degenerate side-only triangle inequality. This does not change the
  separate side-length-based `Triangle` class or promise exact orientation
  predicates or arbitrary extreme-scale arithmetic.
  The height must still be present in the stored `Real` coordinates: if adding
  a translation rounds it away, the resulting collinear vertices have zero
  area and neither winding. Choosing a local coordinate frame can avoid that
  input-rounding loss; increasing an area tolerance cannot recover lost data.
- `TriangleSurface3D` publicly implements `IParametricSurface<3>`, including
  pointer/reference conversion and polymorphic evaluation through its inherited
  `IVectorFunctionNM<2, 3>` interface. Existing parameter bounds, vertex
  reordering and local-coordinate mapping are unchanged.
- `Line3D::Intersection` and both line-line `Dist` overloads use the existing
  `IsParallel` direction tolerance consistently. Very-nearly-parallel lines are
  deliberately treated as a parallel approximation, even if their infinite
  extensions meet far away. For this branch, detailed results fix `t2 = 0`,
  project the second reference point onto the first line, and populate both
  points, parameters and separation, including for coincident lines. This
  representative pair is not a unique minimum for tolerance-parallel directions.
  The distance overload with point outputs retains its `false` return and leaves
  those outputs unchanged for parallel/coincident lines; its distance is populated.
  Nonparallel pairs use a cross-product solve rather than a cancellation-prone
  Gram determinant. Separation is evaluated from the line geometry, not by
  subtracting potentially large closest-point coordinates. Invalid inputs and
  unrepresentable required results throw `GeometryError`; scalar distance does
  not require representable closest-point parameters.
- `SegmentLine2D::ClosestPoint` and `DistanceToPoint` treat only identical
  endpoints as a point-segment. Nonzero small segments retain their projections;
  scaled projection and `hypot` avoid squaring tiny or large lengths.
  Nonfinite coordinates or unrepresentable required coordinate differences and
  distances throw `GeometryError`. Intersections involving a point-segment test
  its actual distance to the other segment against `Precision::GeometryEpsilon`,
  in either operand order. Two point-segments intersect within that same distance
  tolerance. Nondegenerate orientation/collinearity tolerances are unchanged.
- `Tetrahedron` orders its edge lengths as AB, AC, AD, BC, BD, CD.
  `Volume()` validates the current lengths, including changes through mutable accessors.
  Nonfinite or negative lengths, violated face triangle inequalities, nonrealizable
  tetrahedra, and volume overflow throw `GeometryError`. Realizable degenerate shapes
  have zero volume. Validation uses scaled roundoff allowances; only roundoff-sized
  negative Gram determinants are clamped to zero, not arbitrary negative determinants.
  Range tests must keep the cubic volume representable in the active `Real`
  type; representable edge lengths alone do not guarantee a representable volume.
- `Parallelogram::FromDiagonals(p, q, angleBetween)` takes positive diagonal lengths
  and an angle strictly between zero and pi. The returned `DiagonalP()` and
  `DiagonalQ()` match `p` and `q`; its included side angle follows that diagonal
  convention and need not equal `angleBetween`.
  Scaled half-angle formulas and `atan2` avoid cancellation in side/angle
  reconstruction; the diagonal accessors also use half-angle norms rather than
  subtracting nearly equal squares. The stored `Real` angle still limits
  resolution for included angles extremely close to pi.
- A plane-plane intersection returns a line whose reference point satisfies both
  plane equations, not merely a direction parallel to both planes. Parallel planes
  return `false`.
- `Pyramid3D` places the square base in the XY plane through its base center and
  its apex one height above that point. Both constructors generate the same geometry
  up to translation. Its six outward triangular faces agree with analytic volume,
  area, centroid and bounding box.
- `PyramidEquilateral3D` is a **square-based pyramid with all eight edges equal**,
  not a regular tetrahedron. Its height is the base side length divided by sqrt(2).
- Sphere triangles use coherent outward winding, including both pole caps.
  Sphere and torus `Volume()` and `SurfaceArea()` retain their smooth-body analytic
  formulas; their finite meshes approximate those curved bodies.
  Closed triangle/quad mesh centroids accumulate first moments relative to a
  mesh vertex, then restore the translation. Large world coordinates therefore
  do not enter the moment sums; input vertex quantization and ordinary summation
  roundoff still limit accuracy.
- `RectSurface3D` uses a bilinear map over the existing rectangular parameter bounds
  for coplanar, cyclically ordered convex quadrilaterals. Every parameter corner
  maps to the corresponding stored corner, and `getArea()` measures the planar
  quadrilateral rather than an assumed parallelogram. Rectangle and parallelogram
  parameterizations are unchanged. `getCenter()` is the vertex average (the
  parameter-center point), not generally the area centroid. Torus patches retain
  their quad API and face count and use this same mapping and area.
- Point-normal and three-point `Plane3D` construction normalize nonzero directions
  independently of their physical magnitude. Scaled, compensated cross products
  avoid premature product overflow, underflow and cancellation. Nonfinite inputs,
  zero normals/collinear points and unrepresentable plane offsets throw `GeometryError`.
- `RectSurface3D` checks coplanarity relative to the larger of its two adjacent
  edge lengths, using `Precision::Plane3DIsPointOnPlaneTolerance`. This intentionally
  replaces the old absolute-distance constructor check: uniformly scaled copies
  receive the same planarity decision, within floating-point rounding. The public
  `Plane3D::IsPointOnPlane` API retains its **absolute** tolerance.
  Adjacent edge lengths and parameter half-widths must be nonzero and representable.
  `getArea()` retains representable tiny and large areas and throws `GeometryError`
  if the area is zero, underflows to zero or overflows. Centers and adaptive
  rectangular-mesh edge midpoints avoid overflowing coordinate sums.

### Plane query range and coefficient contracts

Plane queries accept finite nonzero normals, including nonunit coefficients from
the intercept constructor or mutable `A()`, `B()`, `C()`, `D()` accessors.
`Normal()` returns those stored coefficients, not a newly normalized vector.
Distances, projections, geometric angle/parallel tests and intersections account
for their magnitude. Line/plane intersection parallel tolerances are angular
(normalized normal/direction dot products).

Membership deliberately retains the **absolute algebraic-residual** contract:
`abs(A*x+B*y+C*z+D) < tolerance`, not a coordinate-relative or normalized-distance
test. Rescaling coefficients can therefore change approximate membership.
The strict inequality is unchanged; zero tolerance always returns false.
Line membership evaluates the residual at parameters zero and one without
constructing a potentially overflowing second point.

Shared exponent-separated, compensated arithmetic avoids premature product,
sum, normal-square and intersection-parameter range loss. Projection also
retains cancellation in its coordinate numerators. These are floating-point
queries, not arbitrary-precision output or exact incidence predicates.
Rounding stored offsets or returned coordinates can still matter at a tight
absolute membership tolerance, especially far from the origin.

As approved for these queries, invalid geometry (zero/nonfinite normals,
nonfinite offsets/points or invalid line directions) raises `GeometryError`.
Explicit tolerances must be finite and nonnegative. A genuinely unrepresentable
numeric output, including nonzero underflow to zero, also raises `GeometryError`.
An unrepresentable intermediate alone does not: projection can succeed when
the point-to-plane distance is too large to return, and an intersection can
succeed when its line parameter is outside `Real` range.
Valid parallel intersections still return false; axis-parallel intercepts
still return positive infinity by definition. Output arguments are unchanged
on false or an exception.

`getPointOnPlane()` now deterministically returns the closest point to the
origin, rather than selecting an arbitrary coordinate-axis intercept.
This can remain finite when every axis intercept is out of range.

This checked example is runnable through `MML_DocsApp --plane-query-range`:

```cpp
const Real largest = std::numeric_limits<Real>::max();
const Plane3D plane(Pnt3Cart(0, 0, 0), Vec3Cart(1, 1, 1));
const Pnt3Cart point(largest, largest, -largest / 2);
const Real distance = plane.DistToPoint(point);
const auto projection = plane.ProjectionToPlane(point);
std::cout << "Scaled plane distance: " << distance / largest << '\n';
std::cout << "Scaled projection: " << projection.X() / largest << ", "
    << projection.Y() / largest << ", " << projection.Z() / largest << '\n';
```

The following construction example is runnable through `MML_DocsApp --plane-rectangle-contracts`:

```cpp
const Real side = Real(1e-8);
const RectSurface3D square(Pnt3Cart(0, 0, 0), Pnt3Cart(side, 0, 0),
    Pnt3Cart(side, side, 0), Pnt3Cart(0, side, 0));
std::cout << "Scaled patch area: " << square.getArea() << '\n'; // 1e-16
```

---

## Examples

### Example 1: 2D Distance Calculation
```cpp
Point2Cartesian a(0, 0);
Point2Cartesian b(3, 4);

Real dist = a.Dist(b);  // 5.0 (Pythagorean triple)

std::cout << "Distance: " << dist << std::endl;
```

### Example 2: Polar to Cartesian Conversion
```cpp
// Point at 45° angle, distance 5 from origin
Point2Polar polar(5.0, Constants::PI / 4);

// Convert to Cartesian
Point2Cartesian cart = polar.TransfToCart();
// cart ≈ (3.536, 3.536)

std::cout << "Cartesian: (" << cart.X() << ", " << cart.Y() << ")" << std::endl;
```

### Example 3: 3D Point Arithmetic
```cpp
Point3Cartesian origin(0, 0, 0);
Point3Cartesian offset(10, 20, 30);

// Translate by offset
Point3Cartesian p1(5, 5, 5);
Point3Cartesian p2 = p1 + offset;  // (15, 25, 35)

// Scale from origin
Point3Cartesian p3 = p1 * 2.0;  // (10, 10, 10)

// Midpoint
Point3Cartesian mid = (p1 + p2) / 2.0;
```

### Example 4: Spherical Coordinates
```cpp
// Point on unit sphere at (θ=90°, φ=45°)
Point3Spherical sph(1.0, Constants::PI/2, Constants::PI/4);

// Convert to Cartesian
Point3Cartesian cart = sph.TransfToCart();
// cart ≈ (0.707, 0.707, 0)

std::cout << "Cartesian: (" << cart.X() << ", " 
          << cart.Y() << ", " << cart.Z() << ")" << std::endl;
```

### Example 5: Cylindrical Coordinates
```cpp
// Point at radius 5 from z-axis, angle 30°, height 10
Point3Cylindrical cyl(5.0, Constants::PI/6, 10.0);

// Convert to Cartesian
Point3Cartesian cart = cyl.TransfToCart();
// cart ≈ (4.330, 2.5, 10)

// Distance from origin
Real dist = Point3Cartesian(0,0,0).Dist(cart);
```

---

## Type Aliases

```cpp
// Common aliases for brevity
typedef Point3Cartesian    Pnt3Cart;
typedef Point3Spherical    Pnt3Sph;
typedef Point3Cylindrical  Pnt3Cyl;
```

---

## Best Practices

1. **Choose coordinate system based on symmetry:**
   - Cartesian: General purpose, no symmetry
   - Polar/Cylindrical: Rotational symmetry around z-axis
   - Spherical: Radial symmetry from origin

2. **Tolerance for comparisons:**
   - Use `IsEqual(other, eps)` for approximate equality
   - Default tolerances in `Defaults` namespace
   - Adjust `eps` based on problem scale

3. **Conversions:**
   - Convert to Cartesian for general arithmetic
   - Use native coordinate system for problem-specific operations

4. **Distance calculations:**
   - Direct `.Dist()` methods available in each coordinate system
   - For polar: uses law of cosines
   - For spherical/cylindrical: converts to Cartesian internally

---

## See Also
- [Geometry_2D_3D.md](Geometry_2D_3D.md) - Lines, planes, triangles, polygons
- [Vectors.md](Vectors.md) - Vector types for directions and displacements
- [Coordinate_transformations.md](../core/Coordinate_transformations.md) - Systematic coordinate transformations

---

## Runnable Examples

| Example | Description | Source |
|---------|-------------|--------|
| 2D Points Demo | Point2Cartesian and Point2Polar operations | [`docs_demo_geometry.cpp`](../../src/docs_demos/docs_demo_geometry.cpp) |
| 3D Points Demo | Point3Cartesian, Spherical, and Cylindrical | [`docs_demo_geometry.cpp`](../../src/docs_demos/docs_demo_geometry.cpp) |

**To run:** Build and execute `MML_DocsApp` target.
