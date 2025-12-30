# MML Base Layer Documentation

## Overview

The Base layer provides foundational mathematical types and structures that serve as building blocks for all higher-level operations in MinimalMathLibrary. These are template-based, header-only implementations designed for efficiency, type safety, and mathematical correctness.

**Files**: `mml/base/*.h`

### Highlights and Quick Navigation
- **Vectors**: [Vectors.md](Vectors.md) - Overview, [Vector.md](Vector.md) - Dynamic, [VectorN.md](VectorN.md) - Fixed-size
- **Matrices**: [Matrices.md](Matrices.md) - Overview, [Matrix.md](Matrix.md) - Dynamic, [MatrixNM.md](MatrixNM.md) - Fixed-size
- **Tensors**: [Tensors.md](Tensors.md)
- **Geometry**: [Geometry.md](Geometry.md), [Geometry_2D_3D.md](Geometry_2D_3D.md)
- **Polynomials**: [Polynoms.md](Polynoms.md)
- **Intervals**: [Intervals.md](Intervals.md)
- **Standard Functions**: [Standard_functions.md](Standard_functions.md)
- **Utilities**: [BaseUtils.md](BaseUtils.md)

## Contents

1. [Linear Algebra Types](#linear-algebra-types)
   - [Vector](#vector)
   - [VectorN](#vectorn)
   - [Matrix](#matrix)
   - [Specialized Matrices](#specialized-matrices)
   - [Tensor](#tensor)
2. [Geometric Types](#geometric-types)
   - [Geometry2D](#geometry2d)
   - [Geometry3D](#geometry3d)
   - [GeometrySpherical](#geometryspherical)
   - [Quaternions](#quaternions)
3. [Function Types](#function-types)
   - [Function Wrappers](#function-wrappers)
   - [InterpolatedFunction](#interpolatedfunction)
   - [StandardFunctions](#standardfunctions)
4. [Mathematical Structures](#mathematical-structures)
   - [Intervals](#intervals)
   - [Polynom](#polynom)
5. [Differential Equations](#differential-equations)
   - [ODESystem](#odesystem)
   - [ODESystemSolution](#odesystemsolution)
6. [Utilities](#utilities)
   - [BaseUtils](#baseutils)
   - [VectorTypes](#vectortypes)
   - [Random](#random)
   - [DiracDeltaFunction](#diracdeltafunction)

---

## Linear Algebra Types

### Vector

**File**: `mml/base/Vector.h`

**Purpose**: Dynamic-size vector with efficient storage and comprehensive operations.

#### Key Features

- Dynamic sizing with move semantics
- Type-generic (works with `double`, `float`, `complex<double>`, etc.)
- STL-compatible iterators
- Exception-safe memory management
- Arithmetic operations (`+`, `-`, `*`, `/`, scalar multiplication)
- Vector operations (dot product, cross product, norms)
- Element-wise operations

#### Core API

```cpp
// Construction
Vector<double> v1;                          // Empty vector
Vector<double> v2(5);                       // Size 5, zero-initialized
Vector<double> v3(5, 3.14);                 // Size 5, all values = 3.14
Vector<double> v4({1.0, 2.0, 3.0});        // Initializer list
Vector<double> v5 = Vector<double>::GetUnitVector(3, 1);  // e₁ = [0,1,0]

// Access
v[i]                    // Element access with bounds checking in debug
v.front(), v.back()     // First and last elements
v.size()                // Vector dimension

// Iterators
for (auto& val : v) { }  // Range-based for loop

// Operations
v1 + v2                 // Vector addition
v1 - v2                 // Vector subtraction
v1 * 2.0                // Scalar multiplication
v.NormL2()              // Euclidean norm: √(Σxᵢ²)
v.NormL1()              // L1 norm: Σ|xᵢ|
v.NormLInf()            // Infinity norm: max|xᵢ|

// Modification
v.push_back(val)        // Add element
v.Resize(newSize)       // Change size
v.Clear()               // Remove all elements

// Note: Vector<T> has no Normalize() method.
// For unit vector: Vector<double> unit = v / v.NormL2();
// For dot/cross products, use Vector3Cartesian (see VectorTypes.h)
```

#### Usage Example

```cpp
Vector<double> position({1.0, 2.0, 3.0});
Vector<double> velocity({0.5, -0.3, 0.1});

// Calculate displacement
Vector<double> displacement = velocity * deltaTime;

// Update position
position = position + displacement;

// Check if moving
double speed = velocity.NormL2();
if (speed > 0.001) {
    Vector<double> direction = velocity / speed;  // Manual normalization
    // direction is now unit vector in movement direction
}
```

---

### VectorN

**File**: `mml/base/VectorN.h`

**Purpose**: Fixed-size vector with compile-time dimension for maximum performance.

#### Key Features

- Compile-time size (template parameter `N`)
- Stack-allocated (no dynamic memory)
- Optimized for small vectors (2D, 3D, 4D)
- All operations inlined
- Swizzling support for 2D/3D/4D
- Cache-friendly memory layout

#### Core API

```cpp
// Construction
VectorN<double, 3> v1;                     // Zero-initialized
VectorN<double, 3> v2{1.0, 2.0, 3.0};     // List initialization
VectorN<double, 3> v3 = VectorN<double, 3>::GetUnitVector(2);  // e₂

// Access
v[i]                    // Element access
v.size()                // Returns N (constexpr)

// Common operations (same interface as Vector)
v1 + v2, v1 - v2, v1 * scalar
v1.NormL2(), v1.NormLInf()
v.Normalized()          // Returns new unit vector

// Swizzling (for N=2,3,4)
v.x(), v.y(), v.z(), v.w()     // Component access

// Note: For dot/cross products, use Vector2Cartesian/Vector3Cartesian
// which provide ScalarProduct() and VectorProduct() methods
```

#### Usage Example

```cpp
// 3D physics calculation
VectorN<double, 3> force{0.0, -9.81, 0.0};  // Gravity
VectorN<double, 3> velocity{10.0, 5.0, 0.0};

double mass = 2.0;
VectorN<double, 3> acceleration = force / mass;

// No heap allocation - all on stack!
```

#### When to Use

- **Use `VectorN`** when dimension is known at compile time (2D/3D graphics, physics)
- **Use `Vector`** when dimension varies at runtime (general algorithms, data processing)

---

### Matrix

**File**: `mml/base/Matrix.h`

**Purpose**: General-purpose dense matrix with comprehensive linear algebra operations.

#### Key Features

- Dynamic sizing (rows × columns)
- Contiguous memory layout for cache efficiency
- Row-major storage
- Exception-safe allocation (RAII)
- Memory limits to prevent OOM crashes
- Move semantics for efficient returns
- Comprehensive arithmetic and algebraic operations

#### Core API

```cpp
// Construction
Matrix<double> A;                         // Empty matrix
Matrix<double> B(m, n);                   // m×n matrix, zero-initialized
Matrix<double> C(m, n, value);            // m×n matrix, all elements = value
Matrix<double> I = Matrix<double>::GetUnitMatrix(n);     // n×n identity

// Special constructors
Matrix<double> D = Matrix<double>::MakeRowVector(vec);   // 1×n from Vector
Matrix<double> E = Matrix<double>::MakeColVector(vec);   // n×1 from Vector

// Access
A(i, j)                 // Element at row i, column j
A[i]                    // Row i as pointer (use sparingly)
A.RowNum(), A.ColNum()  // Dimensions
A.Print(std::cout, 2)   // Print with 2 decimal places

// Matrix arithmetic
A + B, A - B            // Matrix addition/subtraction
A * B                   // Matrix multiplication
A * scalar              // Scalar multiplication
A.GetTranspose()        // Returns transposed matrix
A.GetInverse()          // Returns inverse matrix (throws if singular)

// In-place operations
A.Transpose()           // Transpose in-place (square matrices only)
A.Invert()              // Invert in-place (throws if singular)

// Properties
A.Trace()               // Trace (sum of diagonal)
// Note: Det, Rank, Norm, ConditionNumber are computed via solver classes
// (LUSolver provides determinant, QRSolver provides rank, etc.)

// Note: Matrix decompositions are performed via separate solver classes:
// LUSolver, QRSolver, CholeskySolver, SVDecompositionSolver
// See Core Layer documentation for linear equation solvers

// Row/Column operations
A.SwapRows(i, j)
A.SwapColumns(i, j)
A.GetRow(i), A.GetColumn(j)
A.SetRow(i, rowVec), A.SetColumn(j, colVec)

// Submatrices
A.GetSubmatrix(row_start, row_end, col_start, col_end)
```

#### Usage Example

```cpp
// Linear system: Ax = b
Matrix<double> A(3, 3);
A(0,0)=2; A(0,1)=-1; A(0,2)=0;
A(1,0)=-1; A(1,1)=2; A(1,2)=-1;
A(2,0)=0; A(2,1)=-1; A(2,2)=2;

Vector<double> b({1.0, 0.0, 1.0});

// Solve using LU decomposition (via solver class)
LUSolver<double> solver(A);
Vector<double> x = solver.Solve(b);

// Verify: should get b back
Vector<double> check = A * x;
double error = (check - b).NormL2();  // Should be ~0

// Matrix operations
Matrix<double> AtA = A.GetTranspose() * A;   // Normal equations
// Note: Condition number available via SVDecompositionSolver
```

#### Memory Safety

```cpp
// Protected against overflow and OOM
try {
    Matrix<double> huge(1000000, 1000000);  // Throws std::bad_alloc
} catch (const std::exception& e) {
    // Gracefully handle allocation failure
}
```

---

### Specialized Matrices

#### MatrixSym

**File**: `mml/base/MatrixSym.h`

**Purpose**: Symmetric matrix with optimized storage (stores only upper triangle).

```cpp
MatrixSym<double> S(n);      // Symmetric n×n matrix
// Uses n(n+1)/2 storage instead of n²
// Automatically maintains symmetry: S(i,j) = S(j,i)
```

**Use Cases**: Covariance matrices, Hessians, quadratic forms

#### MatrixTriDiag

**File**: `mml/base/MatrixTriDiag.h`

**Purpose**: Tridiagonal matrix with O(n) storage and specialized solvers.

```cpp
MatrixTriDiag<double> T(n);
T.SetDiag(i, value);         // Main diagonal
T.SetSubDiag(i, value);      // Subdiagonal
T.SetSuperDiag(i, value);    // Superdiagonal

Vector<double> x = T.SolveTridiag(b);  // O(n) solver
```

**Use Cases**: Finite difference methods, spline interpolation, heat equation

#### MatrixBandDiag

**File**: `mml/base/MatrixBandDiag.h`

**Purpose**: Banded matrix with customizable bandwidth.

```cpp
MatrixBandDiag<double> B(n, lower_bw, upper_bw);
// Stores only the band - massive savings for large sparse systems
```

**Use Cases**: Finite element methods, differential operators on grids

#### Matrix3D

**File**: `mml/base/Matrix3D.h`

**Purpose**: 3×3 matrix optimized for 3D transformations.

```cpp
Matrix3D<double> R = Matrix3D<double>::RotationAroundX(angle);
Matrix3D<double> S = Matrix3D<double>::ScaleMatrix(sx, sy, sz);

VectorN<double, 3> v{1,0,0};
VectorN<double, 3> transformed = R * v;

double det = R.Det3x3();     // Fast 3×3 determinant
```

**Use Cases**: 3D graphics, robotics, physics transformations

#### MatrixNM

**File**: `mml/base/MatrixNM.h`

**Purpose**: Fixed-size matrix (compile-time dimensions).

```cpp
MatrixNM<double, 3, 4> M;    // 3×4 matrix, stack-allocated
// All operations compile-time optimized
// No heap allocation
```

**Use Cases**: Small matrices with known size, performance-critical code

---

### Tensor

**File**: `mml/base/Tensor.h`

**Purpose**: Multi-dimensional array for rank-2 and rank-3 tensors.

#### Key Features

- Support for rank-2 (matrix-like) and rank-3 tensors
- Index notation support
- Tensor contraction operations
- Einstein summation convention

#### Core API

```cpp
// Rank-2 tensor (generalized matrix)
Tensor2<double> T2(n, m);
T2(i, j) = value;

// Rank-3 tensor
Tensor3<double> T3(n, m, k);
T3(i, j, k) = value;

// Tensor operations
T2.Contract(indices);        // Tensor contraction
T2.Transpose(idx1, idx2);    // Generalized transpose
```

#### Usage Example

```cpp
// Stress tensor in continuum mechanics
Tensor2<double> stress(3, 3);
stress(0,0) = sigma_xx;
stress(0,1) = stress(1,0) = tau_xy;  // Symmetry
stress(0,2) = stress(2,0) = tau_xz;
// ... etc

// Von Mises stress calculation
double vonMises = CalculateVonMisesStress(stress);
```

---

## Geometric Types

### Geometry2D

**File**: `mml/base/Geometry2D.h`

**Purpose**: 2D geometric primitives and operations.

#### Components

```cpp
// Points
Point2Cartesian p{x, y};

// Lines
Line2D line{point, direction};

// Circles
Circle2D circle{center, radius};

// Polygons
Polygon2D poly{vertices};

// Operations
bool intersects = LineCircleIntersection(line, circle);
double dist = PointToLineDistance(point, line);
bool inside = PointInPolygon(point, poly);
```

**Use Cases**: Computational geometry, 2D physics, graphics

---

### Geometry3D

**File**: `mml/base/Geometry3D.h`

**Purpose**: 3D geometric primitives and spatial operations.

#### Components

```cpp
// Points and vectors
Point3Cartesian p{x, y, z};
Vector3Cartesian v{vx, vy, vz};

// Lines and rays
Line3D line{point, direction};
Ray3D ray{origin, direction};

// Planes
Plane3D plane{point, normal};

// Spheres and boxes
Sphere3D sphere{center, radius};
Box3D box{min_corner, max_corner};

// Triangles
Triangle3D tri{p1, p2, p3};

// Operations
bool hit = RaySphereIntersection(ray, sphere, hitPoint);
double dist = PointToPlaneDistance(point, plane);
VectorN<double, 3> normal = TriangleNormal(tri);
bool inside = PointInBox(point, box);
```

**Use Cases**: 3D graphics, ray tracing, collision detection, CAD

---

### Geometry3DBodies

**File**: `mml/base/Geometry3DBodies.h`

**Purpose**: 3D solid body primitives with geometric properties for physics simulations, CAD, and collision detection.

#### Body Types

| Body | Description | Formula References |
|------|-------------|-------------------|
| `Sphere3D` | Sphere with radius R | V = (4/3)πR³, SA = 4πR² |
| `Cube3D` | Axis-aligned cube with side a | V = a³, SA = 6a² |
| `Cylinder3D` | Right circular cylinder | V = πR²H, SA = 2πR(R+H) |
| `Torus3D` | Torus with major radius R, minor radius r | V = 2π²Rr², SA = 4π²Rr |
| `Pyramid3D` | Square-base pyramid | V = (1/3)a²h |
| `CubeWithTriangles3D` | Triangulated cube (12 triangles) | Same as Cube3D |
| `PyramidEquilateral3D` | Regular tetrahedron (equilateral) | h = a/√3 constraint |

#### IBody Interface

All bodies implement the `IBody` interface:

```cpp
class IBody {
public:
    virtual Real Volume() const = 0;
    virtual Real SurfaceArea() const = 0;
    virtual Pnt3Cart GetCenter() const = 0;
    virtual BoundingBox3D GetBoundingBox() const = 0;
    virtual BoundingSphere3D GetBoundingSphere() const = 0;
    virtual bool IsInside(const Pnt3Cart& pnt) const = 0;
    virtual std::string ToString() const = 0;
};
```

#### Bounding Volumes

Two bounding volume types for spatial queries and collision detection:

```cpp
// Axis-Aligned Bounding Box (AABB)
BoundingBox3D bbox(minPoint, maxPoint);
bool hit = bbox.Contains(point);
bool overlap = bbox.Intersects(otherBox);
Real vol = bbox.Volume();

// Bounding Sphere (faster intersection tests)
BoundingSphere3D bsphere(center, radius);
bool inside = bsphere.Contains(point);
bool collide = bsphere.Intersects(otherSphere);
```

#### Usage Examples

```cpp
// Create bodies
Sphere3D sphere(5.0);                          // Radius 5, centered at origin
Sphere3D sphere2(3.0, Pnt3Cart(10, 20, 30));   // Radius 3, at specific center

Cube3D cube(10.0);                             // Side 10, centered at origin
Cylinder3D cylinder(2.0, 8.0);                 // Radius 2, height 8
Torus3D torus(5.0, 1.0);                       // Major radius 5, minor radius 1
Pyramid3D pyramid(4.0, 6.0);                   // Base side 4, height 6

// Query geometric properties
Real v = sphere.Volume();           // (4/3)π·125 ≈ 523.599
Real sa = sphere.SurfaceArea();     // 4π·25 ≈ 314.159
Pnt3Cart c = sphere.GetCenter();    // (0, 0, 0)

// Bounding volumes for collision detection
BoundingBox3D bbox = cube.GetBoundingBox();
BoundingSphere3D bsphere = cube.GetBoundingSphere();

// Point-in-body test
bool inside = sphere.IsInside(Pnt3Cart(1, 1, 1));  // true

// String output for debugging
std::cout << sphere.ToString();
// Sphere3D{Center=(0, 0, 0), Radius=5, Volume=523.599, SurfaceArea=314.159}
```

#### Surface-Based Bodies

For surface integration and flux calculations, use triangulated or rectangular surface bodies:

```cpp
// Bodies with triangular surfaces (for general meshes)
class BodyWithTriangleSurfaces : public IBody {
    int GetSurfaceCount() const;
    const Triangle3D& GetSurface(int index) const;
};

// Bodies with rectangular surfaces (for boxes, etc.)
class BodyWithRectSurfaces : public IBody {
    int GetSurfaceCount() const;
    const RectSurface3D& GetSurface(int index) const;
};

// Example: Surface integration over a cube
CubeWithTriangles3D meshCube(1.0);
for (int i = 0; i < meshCube.GetSurfaceCount(); i++) {
    const Triangle3D& tri = meshCube.GetSurface(i);
    // Process each triangle surface
}
```

**Use Cases**: Physics simulations, collision detection, CAD modeling, surface/volume integration, bounding volume hierarchies

---

### GeometrySpherical

**File**: `mml/base/GeometrySpherical.h`

**Purpose**: Spherical coordinate geometry for astronomy and physics.

#### Key Features

- Spherical coordinates (r, θ, φ)
- Great circle distances
- Spherical trigonometry
- Astronomical calculations

```cpp
Point3Spherical p{r, theta, phi};

// Convert to/from Cartesian
Point3Cartesian cart = SphericalToCartesian(p);
Point3Spherical sph = CartesianToSpherical(cart);

// Great circle distance
double angularDist = GreatCircleDistance(p1, p2);

// Spherical law of cosines
double sideC = SphericalLawOfCosines(sideA, sideB, angleC);
```

**Use Cases**: Astronomy, geodesy, global positioning, spherical geometry

---

### Quaternions

**File**: `mml/base/Quaternions.h`

**Purpose**: Quaternion algebra for 3D rotations (superior to Euler angles).

#### Key Features

- Gimbal lock-free rotations
- Smooth interpolation (SLERP)
- Efficient composition
- Conversion to/from rotation matrices

#### Core API

```cpp
// Construction
Quaternion q1;                              // Identity (no rotation)
Quaternion q2{w, x, y, z};                 // Direct components
Quaternion q3 = Quaternion::FromAxisAngle(axis, angle);
Quaternion q4 = Quaternion::FromEuler(roll, pitch, yaw);
Quaternion q5 = Quaternion::FromRotationMatrix(R);

// Operations
q1 * q2                     // Quaternion multiplication (rotation composition)
q1.Conjugate()              // Conjugate (inverse for unit quaternions)
q1.Inverse()                // Full inverse
q1.Normalize()              // Make unit quaternion

// Rotation application
VectorN<double, 3> v{1, 0, 0};
VectorN<double, 3> rotated = q.Rotate(v);

// Interpolation
Quaternion interp = Quaternion::Slerp(q1, q2, t);  // t ∈ [0,1]

// Conversions
Matrix3D<double> R = q.ToRotationMatrix();
Vector<double> euler = q.ToEulerAngles();

// Properties
q.Norm()                    // Quaternion magnitude
q.Angle()                   // Rotation angle
VectorN<double, 3> axis = q.Axis();  // Rotation axis
```

#### Usage Example

```cpp
// Camera rotation system
Quaternion cameraOrientation;  // Start at identity

// Rotate around Y axis (yaw)
Quaternion yaw = Quaternion::FromAxisAngle({0,1,0}, yawAngle);
cameraOrientation = yaw * cameraOrientation;

// Rotate around local X axis (pitch)
VectorN<double, 3> localX = cameraOrientation.Rotate({1,0,0});
Quaternion pitch = Quaternion::FromAxisAngle(localX, pitchAngle);
cameraOrientation = pitch * cameraOrientation;

// Apply to view direction
VectorN<double, 3> viewDir = cameraOrientation.Rotate({0,0,-1});

// Smooth rotation interpolation for animation
Quaternion target = Quaternion::FromEuler(targetRoll, targetPitch, targetYaw);
for (double t = 0.0; t <= 1.0; t += 0.01) {
    Quaternion current = Quaternion::Slerp(cameraOrientation, target, t);
    RenderFrame(current);
}
```

**Advantages over Euler Angles**:
- No gimbal lock
- Smooth interpolation
- Efficient composition
- Numerically stable

---

## Function Types

### Function Wrappers

**File**: `mml/base/Function.h`

**Purpose**: Type-safe wrappers for mathematical functions with uniform interface.

#### Function Types

##### RealFunction (ℝ → ℝ)

```cpp
// From function pointer
double myFunc(const double x) { return x*x; }
RealFunction f(myFunc);
double y = f(2.0);  // y = 4.0

// From lambda
RealFunctionFromStdFunc g([](const double x) { return std::sin(x); });
double z = g(PI/2);  // z = 1.0
```

##### ScalarFunction (ℝⁿ → ℝ)

```cpp
// Multivariate function
double distance(const VectorN<double, 3>& v) {
    return v.NormL2();
}
ScalarFunction<3> f(distance);
double d = f({3.0, 4.0, 0.0});  // d = 5.0
```

##### VectorFunction (ℝⁿ → ℝⁿ)

```cpp
// Vector field
VectorN<double, 2> rotateField(const VectorN<double, 2>& v) {
    return {-v[1], v[0]};  // 90° rotation
}
VectorFunction<2> F(rotateField);
auto result = F({1.0, 0.0});  // result = {0, 1}
```

##### VectorFunctionNM (ℝⁿ → ℝᵐ)

```cpp
// Transform from R³ to R²
VectorN<double, 2> project(const VectorN<double, 3>& v) {
    return {v[0], v[1]};  // Project to xy-plane
}
VectorFunctionNM<3, 2> proj(project);
```

##### ParametricCurve

```cpp
template<int N>
class ParametricCurve : public IParametricCurve<N> {
    VectorN<Real, N> (*_func)(Real t);
    Real _minT, _maxT;
public:
    ParametricCurve(VectorN<Real, N> (*func)(Real));  // Unbounded
    ParametricCurve(Real minT, Real maxT, VectorN<Real, N> (*func)(Real));
    
    VectorN<Real, N> operator()(Real t) const;
    Real getMinT() const;
    Real getMaxT() const;
};

// Also available: ParametricCurveFromStdFunc (uses std::function)

// Example: Circle
VectorN<Real, 2> circle(Real t) {
    return {std::cos(t), std::sin(t)};
}
ParametricCurve<2> C(0.0, 2*Constants::PI, circle);
auto point = C(Constants::PI/4);  // Point on curve

// Note: Derivatives and arc length computed via Core/Algorithms layers
```

##### ParametricSurface

```cpp
// For non-rectangular domains (y bounds depend on x)
template<int N>
class ParametricSurface : public IParametricSurface<N> {
public:
    ParametricSurface(VectorN<Real, N>(*func)(Real x, Real y),
                      Real minX, Real maxX,
                      Real(*y1)(Real), Real(*y2)(Real));
    
    VectorN<Real, N> operator()(Real x, Real y) const;
    Real getMinU() const;
    Real getMaxU() const;
    Real getMinW(Real x) const;  // Lower y bound at x
    Real getMaxW(Real x) const;  // Upper y bound at x
};

// For rectangular domains
template<int N>
class ParametricSurfaceRect : public IParametricSurfaceRect<N> {
public:
    ParametricSurfaceRect(VectorN<Real, N>(*func)(Real u, Real w),
                          Real minU, Real maxU, Real minW, Real maxW);
    
    VectorN<Real, N> operator()(Real u, Real w) const;
};

// Example: Sphere
VectorN<Real, 3> sphere(Real u, Real v) {
    return {std::cos(u)*std::sin(v), 
            std::sin(u)*std::sin(v), 
            std::cos(v)};
}
ParametricSurfaceRect<3> S(sphere, 0, 2*Constants::PI, 0, Constants::PI);
```

---

### Interpolated Functions

**File**: `mml/base/InterpolatedFunction.h`

**Purpose**: Create continuous functions from discrete data points.

#### Key Features

- Linear interpolation (`LinearInterpRealFunc`)
- Polynomial interpolation (`PolynomInterpRealFunc`)
- Spline interpolation (`SplineInterpRealFunc`)
- Rational interpolation (`RationalInterpRealFunc`)
- Base class `RealFunctionInterpolated` for common interface

```cpp
// Data points
Vector<Real> x_vals({0.0, 1.0, 2.0, 3.0});
Vector<Real> y_vals({0.0, 1.0, 4.0, 9.0});

// Linear interpolation
LinearInterpRealFunc linear(x_vals, y_vals);
double y1 = linear(0.5);   // Linear interpolation

// Polynomial interpolation (uses all points)
PolynomInterpRealFunc poly(x_vals, y_vals);
double y2 = poly(1.5);     // Polynomial fit

// Spline interpolation (smooth cubic spline)
SplineInterpRealFunc spline(x_vals, y_vals);
double y3 = spline(2.5);   // Smooth interpolation

// Rational interpolation (for functions with poles)
RationalInterpRealFunc rational(x_vals, y_vals);
double y4 = rational(0.8); // Rational fit
```

**Use Cases**: Data fitting, experimental data analysis, lookup tables

---

### Standard Functions

**File**: `mml/base/StandardFunctions.h`

**Purpose**: Pre-defined common mathematical functions.

#### Available Functions

```cpp
namespace Functions {
    // Basic functions (templated - work with Real, Complex, float, etc.)
    Sin<T>, Cos<T>, Tan<T>, Sec<T>, Csc<T>, Ctg<T>
    Exp<T>, Log<T>, Log10<T>, Sqrt<T>, Pow<T>
    Sinh<T>, Cosh<T>, Tanh<T>, Sech<T>, Csch<T>, Ctgh<T>
    Asin<T>, Acos<T>, Atan<T>
    Asinh<T>, Acosh<T>, Atanh<T>
    
    // Special functions (Real only)
    Erf(x)           // Error function
    Erfc(x)          // Complementary error function
    TGamma(x)        // Gamma function
    LGamma(x)        // Log-gamma function
    RiemannZeta(x)   // Riemann zeta function
    Beta(x, y)       // Beta function
    Expint(x)        // Exponential integral
    
    // Orthogonal polynomials
    Hermite(n, x)    // Hermite polynomials
    Legendre(n, x)   // Legendre polynomials
    Laguerre(n, x)   // Laguerre polynomials
    AssocLegendre(l, m, x)   // Associated Legendre
    AssocLaguerre(n, m, x)   // Associated Laguerre
    
    // Bessel functions
    CylBesselJ(n, x) // Bessel J (first kind)
    CylBesselY(n, x) // Bessel Y (second kind)
    CylBesselI(n, x) // Modified Bessel I
    CylBesselK(n, x) // Modified Bessel K
    SphBessel(n, x)  // Spherical Bessel
    SphNeumann(n, x) // Spherical Neumann
    
    // Elliptic integrals
    Ellint_1(k, phi)        // Incomplete elliptic integral (1st kind)
    Ellint_2(k, phi)        // Incomplete elliptic integral (2nd kind)
    Comp_ellint_1(k)        // Complete elliptic integral (1st kind)
    Comp_ellint_2(k)        // Complete elliptic integral (2nd kind)
}
```

**Usage**:
```cpp
Real y = Functions::Sin(Constants::PI / 2);  // y = 1.0
Real g = Functions::TGamma(5.0);             // g = 24.0 (4!)
Real j0 = Functions::CylBesselJ(0, 1.0);     // Bessel J₀(1)
```

---

## Mathematical Structures

### Intervals

**File**: `mml/base/Intervals.h`

**Purpose**: Interval arithmetic for representing and manipulating ranges.

**See detailed documentation**: [Intervals.md](Intervals.md)

#### Quick Reference

```cpp
// Simple intervals
ClosedInterval<double> a(0.0, 10.0);      // [0, 10]
OpenInterval<double> b(0.0, 10.0);        // (0, 10)
RightRay<double> c(5.0);                  // [5, +∞)

// Compound intervals (unions)
Interval<double> compound;
compound.addInterval(make_unique<ClosedInterval<double>>(0,5));
compound.addInterval(make_unique<ClosedInterval<double>>(10,15));
// Represents [0,5] ∪ [10,15]

// Set operations
auto intersection = Interval<double>::Intersection(a, b);
auto difference = Interval<double>::Difference(a, b);
auto complement = Interval<double>::Complement(a);

// Queries
bool contains = a.contains(5.0);
double length = a.getLength();

// Generate uniform points
std::vector<double> points;
a.GetEquidistantCovering(10, points);
```

---

### Polynom

**File**: `mml/base/Polynom.h`

**Purpose**: Polynomial representation and operations.

#### Key Features

- Arbitrary degree polynomials
- Arithmetic operations
- Root finding
- Evaluation (Horner's method)
- Differentiation and integration

#### Core API

```cpp
// Construction
Polynom p({1.0, -3.0, 2.0});  // p(x) = x² - 3x + 2

// Evaluation
double y = p(2.0);             // p(2) = 0

// Arithmetic
Polynom sum = p1 + p2;
Polynom product = p1 * p2;
Polynom scaled = p * 3.0;

// Calculus
Polynom derivative = p.Derive();
Polynom integral = p.Integrate();  // + C (constant term = 0)

// Note: Root finding is in Algorithms layer
// See RootFindingPolynoms.h for Laguerre, Durand-Kerner methods

// Properties
int degree = p.GetDegree();
CoefT coeff = p[i];         // Get coefficient at index i
```

#### Usage Example

```cpp
// Quadratic: x² - 5x + 6 = (x-2)(x-3)
Polynom quadratic({6.0, -5.0, 1.0});  // Coefficients: [c₀, c₁, c₂]

// Derivative: 2x - 5
Polynom deriv = quadratic.Derive();

// Evaluate at x = 2.5 (critical point where deriv = 0)
double minValue = quadratic(2.5);  // -0.25

// Note: For root finding, use Algorithms layer:
// auto roots = PolynomLaguerreAllRoots(quadratic);
```

---

## Differential Equations

### ODESystem

**File**: `mml/base/ODESystem.h`

**Purpose**: Representation of systems of Ordinary Differential Equations.

#### Key Features

- Generic ODE system: **dy/dt = f(t, y)**
- Arbitrary system size
- State-space representation
- Compatible with various solvers

#### Structure

```cpp
template<typename Real>
class ODESystem {
    int _dim;                          // System dimension
    Vector<Real> (*_func)(Real, const Vector<Real>&);  // RHS function
    
public:
    ODESystem(int dim, Vector<Real> (*func)(Real, const Vector<Real>&))
        : _dim(dim), _func(func) {}
    
    Vector<Real> operator()(Real t, const Vector<Real>& y) const {
        return _func(t, y);
    }
    
    int Dim() const { return _dim; }
};
```

#### Usage Example

```cpp
// Simple harmonic oscillator: ẍ + ω²x = 0
// Convert to system: y = [x, ẋ]ᵀ, dy/dt = [ẋ, -ω²x]ᵀ

Vector<double> SHO_RHS(double t, const Vector<double>& y) {
    double omega = 1.0;
    Vector<double> dydt(2);
    dydt[0] = y[1];                    // dx/dt = v
    dydt[1] = -omega*omega * y[0];     // dv/dt = -ω²x
    return dydt;
}

ODESystem<double> oscillator(2, SHO_RHS);

// Can now be passed to any ODE solver
```

**Common Patterns**:

```cpp
// Predator-Prey (Lotka-Volterra)
Vector<double> PredatorPrey(double t, const Vector<double>& y) {
    double alpha = 0.1, beta = 0.02, gamma = 0.3, delta = 0.01;
    Vector<double> dydt(2);
    dydt[0] = alpha*y[0] - beta*y[0]*y[1];      // Prey
    dydt[1] = delta*y[0]*y[1] - gamma*y[1];     // Predators
    return dydt;
}

// N-body problem
Vector<double> NBody(double t, const Vector<double>& y) {
    // y = [x₁,y₁,z₁,vx₁,vy₁,vz₁, x₂,y₂,z₂,vx₂,vy₂,vz₂, ...]
    int N = y.size() / 6;
    Vector<double> dydt(6*N);
    // Compute gravitational forces...
    return dydt;
}
```

---

### ODESystemSolution

**File**: `mml/base/ODESystemSolution.h`

**Purpose**: Container for ODE solution trajectories.

#### Key Features

- Stores time series: t₀, t₁, ..., tₙ
- Stores state vectors: y(t₀), y(t₁), ..., y(tₙ)
- Interpolation between steps
- Export to files
- Visualization support

```cpp
class ODESystemSolution {
    std::vector<double> _times;         // Time points
    std::vector<Vector<double>> _states; // State at each time
    
public:
    void AddPoint(double t, const Vector<double>& y);
    Vector<double> Interpolate(double t) const;
    
    void ExportToFile(const std::string& filename) const;
    void Print(std::ostream& out) const;
    
    int NumPoints() const;
    double TimeRange() const;
};
```

#### Usage Example

```cpp
// After solving ODE
ODESystemSolution sol = solver.Solve(system, y0, t0, tf);

// Access solution
for (int i = 0; i < sol.NumPoints(); ++i) {
    double t = sol.GetTime(i);
    Vector<double> y = sol.GetState(i);
    std::cout << "t=" << t << ", x=" << y[0] << ", v=" << y[1] << "\n";
}

// Interpolate between steps
Vector<double> y_mid = sol.Interpolate(2.5);

// Export for plotting
sol.ExportToFile("trajectory.dat");
```

---

## Utilities

### BaseUtils

**File**: `mml/base/BaseUtils.h`

**Purpose**: Utility functions for common mathematical operations.

```cpp
namespace Utils {
    // Levi-Civita symbol
    int LeviCivita(int i, int j, int k);          // 3D permutation symbol
    int LeviCivita(int i, int j, int k, int l);   // 4D permutation symbol
    
    // Angle conversions
    Real DegToRad(Real angleDeg);
    Real RadToDeg(Real angleRad);
    void AngleDegToExplicit(Real angle, Real& deg, Real& min, Real& sec);
    Real ExplicitToAngleDeg(Real deg, Real min, Real sec);
    Real AngleTo2PiRange(Real rad);    // Normalize to [0, 2π)
    Real AngleToPiPiRange(Real rad);   // Normalize to [-π, π)
    
    // Complex number comparison
    bool AreEqual(const Complex& a, const Complex& b, double eps);
    bool AreEqualAbs(const Complex& a, const Complex& b, double eps);
    
    // Vector comparison and operations
    bool AreEqual(const Vector<Real>& a, const Vector<Real>& b, Real eps);
    Real ScalarProduct(const Vector<Real>& a, const Vector<Real>& b);
    Real VectorsAngle(const Vector<Real>& a, const Vector<Real>& b);
    Vector<Real> VectorProjectionParallelTo(const Vector<Real>& orig, const Vector<Real>& b);
    Vector<Real> VectorProjectionPerpendicularTo(const Vector<Real>& orig, const Vector<Real>& b);
    
    // Matrix utilities
    Matrix<Real> ColumnMatrixFromVector(const Vector<Real>& v);
    Matrix<Real> RowMatrixFromVector(const Vector<Real>& v);
}
```

---

### VectorTypes

**File**: `mml/base/VectorTypes.h`

**Purpose**: Specialized 2D and 3D vector/point classes with geometric operations.

```cpp
// 2D Cartesian vector class (extends VectorN<Real, 2>)
class Vector2Cartesian : public VectorN<Real, 2> {
public:
    Vector2Cartesian(Real x, Real y);
    Real X() const;  Real& X();
    Real Y() const;  Real& Y();
    // + all VectorN operations (dot product, norm, etc.)
};

// 3D Cartesian vector class (extends VectorN<Real, 3>)
class Vector3Cartesian : public VectorN<Real, 3> {
public:
    Vector3Cartesian(Real x, Real y, Real z);
    Real X() const;  Real& X();
    Real Y() const;  Real& Y();
    Real Z() const;  Real& Z();
    // + all VectorN operations
};

// Point-Vector interoperability
Point2Cartesian + Vector2Cartesian -> Point2Cartesian
Point3Cartesian + Vector3Cartesian -> Point3Cartesian
```

---

### Random

**File**: `mml/base/Random.h`

**Purpose**: Random number generation for simulations and testing.

```cpp
class Random {
public:
    // Uniform distributions
    static Real UniformReal(Real min, Real max);
    static int UniformInt(int min, int max);
    
    // Random direction vectors (uniform on sphere)
    static Real UniformVecDirection2(Real& vx, Real& vy, Real abs);
    static Real UniformVecDirection3(Real& vx, Real& vy, Real& vz, Real abs);
};
```

#### Usage

```cpp
// Monte Carlo simulation for π
int numSamples = 10000;
int hits = 0;
for (int i = 0; i < numSamples; ++i) {
    Real x = Random::UniformReal(-1.0, 1.0);
    Real y = Random::UniformReal(-1.0, 1.0);
    if (x*x + y*y <= 1.0) hits++;
}
Real pi_estimate = 4.0 * hits / numSamples;

// Random 3D direction with unit length
Real vx, vy, vz;
Random::UniformVecDirection3(vx, vy, vz, 1.0);
```

---

### Dirac Delta Function Approximations

**File**: `mml/base/DiracDeltaFunction.h`

**Purpose**: Approximations of Dirac delta function for numerical work.

```cpp
// Base class
class DiracFunction : public IRealFunction {
protected:
    int _N;    // Width parameter (larger = narrower, more delta-like)
public:
    DiracFunction(int N);
};

// Approximations:
class DiracStep : public DiracFunction {  // Rectangular pulse
    // δ(x) ≈ N for |x| < 1/(2N), 0 otherwise
};

class DiracExp : public DiracFunction {   // Gaussian approximation
    // δ(x) ≈ N/√(2π) * exp(-x²N²/2)
};

class DiracSqr : public DiracFunction {   // Lorentzian approximation
    // δ(x) ≈ N/(π(1 + N²x²))
};

class DiracSin : public DiracFunction {   // Sinc approximation
    // δ(x) ≈ sin(Nx)/(πx)
};
```

**Usage**:
```cpp
DiracExp delta(100);     // Gaussian approximation with N=100
Real y = delta(0.0);     // Peak value at x=0
Real y2 = delta(0.01);   // Value near center
```

**Use Cases**: Signal processing, Green's functions, impulse responses

---

## Design Principles

### Type Safety

All types use strong typing and templates for compile-time safety:
```cpp
Vector<double> v1;
Vector<int> v2;
// v1 + v2;  // Compile error - type mismatch!
```

### Exception Safety

All classes provide strong exception guarantees:
- Constructors either succeed completely or throw (no partial construction)
- Operations maintain valid state even if exceptions occur
- RAII ensures proper cleanup

### Move Semantics

All types support efficient moves:
```cpp
Matrix<double> CreateLargeMatrix() {
    Matrix<double> M(1000, 1000);
    // ... fill matrix ...
    return M;  // Move, not copy!
}

Matrix<double> result = CreateLargeMatrix();  // No copy!
```

### Memory Efficiency

- `VectorN` and `MatrixNM`: Stack allocation
- `Vector` and `Matrix`: Contiguous heap allocation
- Specialized matrices: Optimal storage (symmetric, banded, etc.)

---

## Performance Notes

### Cache Efficiency

- Matrices use row-major, contiguous layout
- Vectors use contiguous `std::vector` storage
- Iterators allow range-based loops with optimal access patterns

### Compile-Time Optimization

- `VectorN`, `MatrixNM`: All operations can inline
- Template specializations for common types
- `constexpr` where applicable

### Runtime Optimization

- Move semantics reduce copies
- Expression templates (planned) for lazy evaluation
- Parallelization opportunities (ThreadPool integration)

---

## Common Patterns

### Linear System Solving

```cpp
Matrix<double> A = /* coefficient matrix */;
Vector<double> b = /* right-hand side */;

// Use solver class from Core layer
LUSolver<double> solver(A);
Vector<double> x = solver.Solve(b);
```

### Eigenvalue Problems

```cpp
Matrix<double> A = /* symmetric matrix */;
// Use EigenSolver from Algorithms layer
// See docs/algorithms/Eigen_solvers.md
```

### Function Minimization Setup

```cpp
ScalarFunction<2> objective(/* your function */);
// Pass to optimization algorithms (Core layer)
```

### ODE System Setup

```cpp
ODESystem<double> system(n, rhs_function);
// Pass to ODE solvers (Algorithms layer)
```

---

## Error Handling

All base types throw meaningful exceptions:

```cpp
try {
    Matrix<double> A(3, 3);
    A.Invert();  // May throw if singular
} catch (const SingularMatrixError& e) {
    std::cerr << "Matrix is singular: " << e.what() << "\n";
} catch (const MatrixDimensionError& e) {
    std::cerr << "Dimension mismatch: " << e.what() << "\n";
}
```

**Common Exception Types**:
- `VectorDimensionError`
- `MatrixDimensionError`
- `SingularMatrixError`
- `InvalidInterval`
- `PolynomialError`

---

## Testing

All base types have comprehensive test coverage in `tests/base/`:

```bash
# Run all base layer tests
cd build
.\tests\Debug\MML_Tests.exe "[base]"

# Run specific tests
.\tests\Debug\MML_Tests.exe "[vector]"
.\tests\Debug\MML_Tests.exe "[matrix]"
.\tests\Debug\MML_Tests.exe "[intervals]"
```

---

## Integration with Higher Layers

The Base layer provides primitives used by:

- **Core Layer**: Operates on these types (differentiation, integration, transformations)
- **Algorithms Layer**: Uses these for numerical methods (solvers, optimizers)
- **Visualizers**: Exports these types for plotting and analysis

**Example Integration**:
```cpp
// Base: Define function
RealFunction f(std::sin);

// Core: Compute derivative
RealFunction df = Derivative(f);

// Algorithms: Find root
double root = NewtonRaphson(f, df, initial_guess);

// Visualizers: Plot
Visualizer::PlotFunction(f, 0.0, 2*PI, "sin(x)");
```

---

## Summary

The Base layer provides:
- **11 core types** for linear algebra (Vector, Matrix, variants)
- **4 geometry systems** (2D, 3D, Spherical, Quaternions)
- **8 function types** (Real, Scalar, Vector, Parametric, Interpolated)
- **2 mathematical structures** (Intervals, Polynomials)
- **ODE system framework** for differential equations
- **Utility functions** for common operations

All types are:
✅ Template-based (type-generic)  
✅ Exception-safe (RAII)  
✅ Move-optimized (zero-copy returns)  
✅ Cache-efficient (contiguous memory)  
✅ Well-tested (comprehensive test suite)  

**Linear Algebra Solvers**: See [Core Layer - Linear Equation Solvers](../core/Linear_equations_solvers.md)  
**Eigenvalue Solvers**: See [Algorithms Layer - Eigen Solvers](../algorithms/Eigen_solvers.md)

**Next Layer**: [Core Layer Documentation](../core/README_Core.md)

---

*Last Updated: December 27, 2025*  
*MML Version: 1.0*
