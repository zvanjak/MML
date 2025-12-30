# Geometry - 2D & 3D Primitives

**Files**: `mml/base/Geometry2D.h`, `mml/base/Geometry3D.h`

Geometric primitives for analytic geometry: lines, segments, planes, triangles, and polygons.

## Table of Contents
- [2D Primitives](#2d-primitives)
- [3D Primitives](#3d-primitives)
- [Examples](#examples)

---

## 2D Primitives

### Line2D

Infinite line in 2D, parameterized as `p(t) = start + t * direction`.

```cpp
// Construction
Point2Cartesian p(0, 0);
Vector2Cartesian dir(1, 1);
Line2D line1(p, dir);  // Point + direction

Point2Cartesian a(0, 0), b(2, 2);
Line2D line2(a, b);    // Two points

// Access
Point2Cartesian start = line1.StartPoint();
Vector2Cartesian direction = line1.Direction();  // Normalized

// Parametric evaluation
Point2Cartesian p_at_2 = line1(2.0);  // Point at t=2
```

### SegmentLine2D

Line segment with defined start and end points.

```cpp
// Construction
Point2Cartesian a(0, 0), b(4, 3);
SegmentLine2D seg1(a, b);

// From point + direction + parameter
Vector2Cartesian dir(3, 4);
SegmentLine2D seg2(a, dir, 5.0);  // Length 5 along dir

// Access
Point2Cartesian start = seg1.StartPoint();
Point2Cartesian end = seg1.EndPoint();
Vector2Cartesian direction = seg1.Direction();
Real length = seg1.Length();

// Parametric point (t ∈ [0, 1])
Point2Cartesian mid = seg1.PointOnSegment(0.5);  // Midpoint
```

### Triangle2D

2D triangle with three vertices.

```cpp
// Construction
Point2Cartesian p1(0, 0), p2(4, 0), p3(2, 3);
Triangle2D tri(p1, p2, p3);

// Access vertices
Point2Cartesian v1 = tri.Pnt1();
Point2Cartesian v2 = tri.Pnt2();
Point2Cartesian v3 = tri.Pnt3();

// Side lengths
Real a = tri.A();  // Distance p1-p2
Real b = tri.B();  // Distance p2-p3
Real c = tri.C();  // Distance p3-p1

// Properties
Real area = tri.Area();          // Heron's formula
bool right = tri.IsRight();      // Right triangle?
bool isosceles = tri.IsIsosceles();    // Two equal sides?
bool equilateral = tri.IsEquilateral(); // All sides equal?
```

### Polygon2D

Simple polygon as ordered list of vertices.

```cpp
// Construction
std::vector<Point2Cartesian> points = {
    {0, 0}, {4, 0}, {4, 3}, {0, 3}
};
Polygon2D poly1(points);

// Initializer list
Polygon2D poly2({Point2Cartesian(0, 0),
                 Point2Cartesian(1, 0),
                 Point2Cartesian(0, 1)});

// Access
std::vector<Point2Cartesian> vertices = poly1.Points();

// Area (shoelace formula)
Real area = poly1.Area();
```

---

## 3D Primitives

### Line3D

Infinite line in 3D: `p(t) = start + t * direction`.

```cpp
// Construction
Point3Cartesian p(0, 0, 0);
Vector3Cartesian dir(1, 0, 0);
Line3D line1(p, dir);

Point3Cartesian a(0, 0, 0), b(1, 1, 1);
Line3D line2(a, b);  // Direction normalized

// Access
Point3Cartesian start = line1.StartPoint();
Vector3Cartesian direction = line1.Direction();

// Parametric evaluation
Point3Cartesian p_at_5 = line1(5.0);

// Queries
bool parallel = line1.IsParallel(line2, 1e-6);
bool perpendicular = line1.IsPerpendicular(line2, 1e-6);
bool equal = (line1 == line2);

// Point queries
Point3Cartesian test(1, 2, 3);
bool on_line = line1.IsPointOnLine(test, 1e-6);
Real distance = line1.Dist(test);
Point3Cartesian nearest = line1.NearestPointOnLine(test);

// Line-line distance
Real line_dist = line1.Dist(line2);
```

### SegmentLine3D

3D line segment.

```cpp
// Construction
Point3Cartesian a(0, 0, 0), b(1, 1, 1);
SegmentLine3D seg1(a, b);

Vector3Cartesian dir(1, 0, 0);
SegmentLine3D seg2(a, dir, 5.0);  // Length 5

// Access
Point3Cartesian start = seg1.StartPoint();
Point3Cartesian end = seg1.EndPoint();
Real length = seg1.Length();
Vector3Cartesian direction = seg1.Direction();

// Parametric point (t ∈ [0, 1])
Point3Cartesian point = seg1.PointOnSegment(0.25);
```

### Plane3D

Plane in 3D: `Ax + By + Cz + D = 0`.

```cpp
// Construction from point + normal
Point3Cartesian p(0, 0, 0);
Vector3Cartesian normal(0, 0, 1);  // z = 0 plane
Plane3D plane1(p, normal);

// From three points
Point3Cartesian p1(1, 0, 0), p2(0, 1, 0), p3(0, 0, 1);
Plane3D plane2(p1, p2, p3);

// From coefficients (Ax + By + Cz + D = 0)
Plane3D plane3(0, 0, 1, 0);  // z = 0

// Predefined planes
auto xy_plane = Plane3D::GetXYPlane();  // z = 0
auto xz_plane = Plane3D::GetXZPlane();  // y = 0
auto yz_plane = Plane3D::GetYZPlane();  // x = 0

// Access
Real a = plane1.A(), b = plane1.B(), c = plane1.C(), d = plane1.D();
Vector3Cartesian n = plane1.Normal();
Point3Cartesian point_on_plane = plane1.GetPointOnPlane();

// Point queries
Point3Cartesian test(1, 2, 3);
bool on_plane = plane1.IsPointOnPlane(test, 1e-6);
Real distance = plane1.DistToPoint(test);
Point3Cartesian projected = plane1.ProjectionToPlane(test);

// Line-plane operations
Line3D line(Point3Cartesian(0, 0, 1), Vector3Cartesian(0, 0, -1));
bool line_on_plane = plane1.IsLineOnPlane(line);
Real angle = plane1.AngleToLine(line);  // Radians

Point3Cartesian intersection;
bool intersects = plane1.IntersectionWithLine(line, intersection);

// Plane-plane operations
Plane3D other(0, 1, 0, 0);  // y = 0
bool parallel = plane1.IsParallelToPlane(other);
bool perpendicular = plane1.IsPerpendicularToPlane(other);
Real plane_angle = plane1.AngleToPlane(other);
Real plane_dist = plane1.DistToPlane(other);

Line3D intersection_line;
bool planes_intersect = plane1.IntersectionWithPlane(other, intersection_line);

// Coordinate axis intersections
Real seg_x, seg_y, seg_z;
plane1.GetCoordAxisSegments(seg_x, seg_y, seg_z);
```

### Triangle3D

3D triangle.

```cpp
// Construction
Point3Cartesian p1(0, 0, 0), p2(1, 0, 0), p3(0, 1, 0);
Triangle3D tri(p1, p2, p3);

// Access
Point3Cartesian v1 = tri.Pnt1();
Point3Cartesian v2 = tri.Pnt2();
Point3Cartesian v3 = tri.Pnt3();
```

---

## Examples

### Example 1: 2D Line-Line Intersection (Manual)
```cpp
Point2Cartesian p1(0, 0), p2(2, 2);
Point2Cartesian p3(2, 0), p4(0, 2);

Line2D line1(p1, p2);  // y = x
Line2D line2(p3, p4);  // y = -x + 2

// Intersection at (1, 1)
// (Manual calculation - parametric form)
```

### Example 2: 3D Line-Plane Intersection
```cpp
// Plane: z = 0 (XY plane)
Plane3D plane = Plane3D::GetXYPlane();

// Line from (0,0,1) going down
Line3D line(Point3Cartesian(0, 0, 1), 
            Vector3Cartesian(0, 0, -1));

// Find intersection
Point3Cartesian hit;
if (plane.IntersectionWithLine(line, hit)) {
    std::cout << "Hit at: (" << hit.X() << ", " 
              << hit.Y() << ", " << hit.Z() << ")" << std::endl;
    // Output: Hit at: (0, 0, 0)
}
```

### Example 3: Plane from Three Points
```cpp
// Define plane through triangle vertices
Point3Cartesian p1(1, 0, 0);
Point3Cartesian p2(0, 1, 0);
Point3Cartesian p3(0, 0, 1);

Plane3D plane(p1, p2, p3);

// Get plane equation
std::cout << "Plane: " << plane.A() << "x + " 
          << plane.B() << "y + " 
          << plane.C() << "z + " 
          << plane.D() << " = 0" << std::endl;
```

### Example 4: Point-to-Line Distance
```cpp
// Line through origin along x-axis
Line3D line(Point3Cartesian(0, 0, 0), 
            Vector3Cartesian(1, 0, 0));

// Point not on line
Point3Cartesian p(0, 3, 4);

// Distance should be 5 (Pythagorean: sqrt(3² + 4²))
Real dist = line.Dist(p);

// Nearest point on line
Point3Cartesian nearest = line.NearestPointOnLine(p);
// nearest = (0, 0, 0)
```

### Example 5: Triangle Area
```cpp
// Right triangle with legs 3 and 4
Point2Cartesian p1(0, 0);
Point2Cartesian p2(3, 0);
Point2Cartesian p3(0, 4);

Triangle2D tri(p1, p2, p3);

Real area = tri.Area();  // 6.0
bool is_right = tri.IsRight();  // true

std::cout << "Area: " << area << std::endl;
std::cout << "Right triangle: " << (is_right ? "yes" : "no") << std::endl;
```

### Example 6: Polygon Area
```cpp
// Square with side length 2
Polygon2D square({
    Point2Cartesian(0, 0),
    Point2Cartesian(2, 0),
    Point2Cartesian(2, 2),
    Point2Cartesian(0, 2)
});

Real area = square.Area();  // 4.0
std::cout << "Square area: " << area << std::endl;
```

### Example 7: Plane-Plane Intersection
```cpp
// XY plane (z = 0)
Plane3D plane1 = Plane3D::GetXYPlane();

// Plane x = 0 (YZ plane)
Plane3D plane2 = Plane3D::GetYZPlane();

// Intersection should be y-axis
Line3D intersection_line;
if (plane1.IntersectionWithPlane(plane2, intersection_line)) {
    std::cout << "Intersection line start: " 
              << intersection_line.StartPoint() << std::endl;
    std::cout << "Direction: " 
              << intersection_line.Direction() << std::endl;
}
```

---

## Best Practices

1. **Tolerance parameters:**
   - Most comparison methods accept `eps` tolerance
   - Use larger tolerances for larger coordinate scales
   - Default tolerances in `Defaults` namespace

2. **Line directions:**
   - Line and segment directions are auto-normalized
   - Parameter `t` has consistent meaning across similar types
   - For segments: `t ∈ [0, 1]`

3. **Plane equations:**
   - Stored as `Ax + By + Cz + D = 0`
   - Normal vector: `(A, B, C)` (may not be normalized)
   - Use `.Normal()` for normalized normal vector

4. **Intersection tests:**
   - Return `bool` indicating success
   - Output parameters (`out_*`) contain results
   - Check return value before using output

5. **Performance:**
   - Fixed-size classes (no dynamic allocation)
   - Inline methods for tight loops
   - Consider VectorN types for batch operations

---

## See Also
- [Geometry.md](Geometry.md) - Point types and coordinate systems
- [Vectors.md](Vectors.md) - Vector operations (cross product, scalar product)
- [Coordinate_transformations.md](../core/Coordinate_transformations.md) - Transform geometry between frames

---

## Runnable Examples

| Example | Description | Source |
|---------|-------------|--------|
| 2D Lines Demo | Line2D and SegmentLine2D operations | [`docs_demo_geometry_2d_3d.cpp`](../../src/docs_demos/docs_demo_geometry_2d_3d.cpp) |
| 2D Triangle Demo | Triangle2D area, type tests | [`docs_demo_geometry_2d_3d.cpp`](../../src/docs_demos/docs_demo_geometry_2d_3d.cpp) |
| 2D Polygon Demo | Polygon2D area, perimeter | [`docs_demo_geometry_2d_3d.cpp`](../../src/docs_demos/docs_demo_geometry_2d_3d.cpp) |
| 3D Lines Demo | Line3D, SegmentLine3D, distances | [`docs_demo_geometry_2d_3d.cpp`](../../src/docs_demos/docs_demo_geometry_2d_3d.cpp) |
| 3D Plane Demo | Plane3D, intersections | [`docs_demo_geometry_2d_3d.cpp`](../../src/docs_demos/docs_demo_geometry_2d_3d.cpp) |
| 3D Triangle Demo | Triangle3D area, centroid | [`docs_demo_geometry_2d_3d.cpp`](../../src/docs_demos/docs_demo_geometry_2d_3d.cpp) |

**To run:** Build and execute `MML_DocsApp` target.


