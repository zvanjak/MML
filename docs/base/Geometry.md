# Geometry - Basic Point Types

**Files**: `mml/base/Geometry.h`, `mml/base/GeometrySpherical.h`

Lightweight point types for 2D and 3D coordinate systems.

## Table of Contents
- [2D Points](#2d-points)
- [3D Points](#3d-points)
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

