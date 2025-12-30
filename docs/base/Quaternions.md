# Quaternions

The `Quaternion` class provides quaternion arithmetic and 3D rotation operations, offering a singularity-free alternative to Euler angles for representing orientations.

**Header:** `base/Quaternions.h`

## Overview

Quaternions extend complex numbers to 4D and are ideal for:
- **3D rotations** without gimbal lock
- **Smooth interpolation** between orientations (SLERP)
- **Compact representation** (4 values vs 9 for rotation matrix)
- **Efficient composition** of multiple rotations

## Mathematical Background

A quaternion is represented as:
```
q = w + xi + yj + zk
```
where:
- `w` is the scalar (real) part
- `(x, y, z)` is the vector (imaginary) part
- `i, j, k` are imaginary units satisfying: `i² = j² = k² = ijk = -1`

**Unit quaternions** (||q|| = 1) represent rotations. A rotation by angle θ around unit axis **u** is:
```
q = [cos(θ/2), sin(θ/2)·u]
```

## Constructors

```cpp
// Default: Identity quaternion (no rotation)
Quaternion q1;                        // [1, 0, 0, 0]

// From components: w + xi + yj + zk
Quaternion q2(0.707, 0.707, 0, 0);    // 90° rotation around X

// From scalar and vector parts
Vec3Cart axis(0, 0, 1);
Quaternion q3(0.707, axis * 0.707);   // 90° around Z

// Pure imaginary (for vector extension)
Quaternion q4(Vec3Cart(1, 2, 3));     // [0, 1, 2, 3]
```

## Factory Methods

### FromAxisAngle
```cpp
// Rotation of 90° (π/2 radians) around Z-axis
Vec3Cart zAxis(0, 0, 1);
auto q = Quaternion::FromAxisAngle(zAxis, Constants::PI / 2);
```

### FromEulerZYX (Yaw-Pitch-Roll)
```cpp
// Aerospace convention: Yaw (Z), Pitch (Y), Roll (X)
Real yaw = 0.5;    // radians
Real pitch = 0.2;
Real roll = 0.1;
auto q = Quaternion::FromEulerZYX(yaw, pitch, roll);
```

### FromEulerXYZ
```cpp
// Alternative convention: Roll (X), Pitch (Y), Yaw (Z)
auto q = Quaternion::FromEulerXYZ(roll, pitch, yaw);
```

### FromRotationMatrix
```cpp
// Convert from 3×3 rotation matrix (uses Shepperd's method)
MatrixNM<Real, 3, 3> rotMat = /* ... */;
auto q = Quaternion::FromRotationMatrix(rotMat);
```

### Identity
```cpp
auto qId = Quaternion::Identity();  // No rotation
```

## Accessors

```cpp
Quaternion q(0.5, 0.5, 0.5, 0.5);

// Component access
Real w = q.w();        // Scalar part: 0.5
Real x = q.x();        // i component: 0.5
Real y = q.y();        // j component: 0.5
Real z = q.z();        // k component: 0.5

// By index: [0]=w, [1]=x, [2]=y, [3]=z
Real w2 = q[0];

// Get scalar and vector parts
Real scalar = q.Scalar();        // 0.5
Vec3Cart vec = q.Vector();       // (0.5, 0.5, 0.5)
```

## Arithmetic Operations

```cpp
Quaternion p(1, 2, 3, 4);
Quaternion q(5, 6, 7, 8);

// Addition / Subtraction
auto sum = p + q;
auto diff = p - q;

// Scalar multiplication / division
auto scaled = p * 2.0;
auto half = p / 2.0;

// Quaternion multiplication (Hamilton product)
// NOTE: Non-commutative! p*q ≠ q*p in general
auto product = p * q;

// In-place operations
p += q;
p *= 2.0;
p *= q;
```

## Quaternion-Specific Operations

### Conjugate and Inverse
```cpp
Quaternion q(0.5, 0.5, 0.5, 0.5);

// Conjugate: q* = w - xi - yj - zk
auto qConj = q.Conjugate();      // [0.5, -0.5, -0.5, -0.5]

// Inverse: q^(-1) = q* / ||q||²
auto qInv = q.Inverse();

// For unit quaternions: Inverse() == Conjugate()
```

### Norm and Normalization
```cpp
Quaternion q(1, 2, 3, 4);

Real normSq = q.NormSquared();   // w² + x² + y² + z² = 30
Real norm = q.Norm();            // √30 ≈ 5.477

// Normalize in place
q.Normalize();                   // Now ||q|| = 1

// Or get normalized copy
Quaternion qNorm = q.Normalized();

// Check if unit quaternion
bool isUnit = q.IsUnit();        // true if ||q|| ≈ 1
bool isId = q.IsIdentity();      // true if q ≈ [1,0,0,0]
```

### Dot Product
```cpp
Quaternion p(1, 0, 0, 0);
Quaternion q(0.707, 0.707, 0, 0);

Real dot = p.Dot(q);  // Useful for interpolation
// dot ≈ 0.707 (cos of half the angle between rotations)
```

## Rotation Operations

### Rotating Vectors
```cpp
// Create rotation: 90° around Z-axis
auto q = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI / 2);

// Rotate a vector
Vec3Cart v(1, 0, 0);             // X-axis
Vec3Cart rotated = q.Rotate(v);  // → (0, 1, 0) Y-axis
```

### Extracting Rotation Parameters
```cpp
// Create a rotation
auto q = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), Constants::PI / 3);

// Get axis and angle back
Vec3Cart axis = q.GetRotationAxis();   // (1, 0, 0)
Real angle = q.GetRotationAngle();     // π/3 ≈ 1.047

// Or both at once
Vec3Cart ax;
Real ang;
q.ToAxisAngle(ax, ang);
```

### Converting to Euler Angles
```cpp
auto q = Quaternion::FromEulerZYX(0.5, 0.2, 0.1);

// Convert back to Euler angles [yaw, pitch, roll]
Vec3Cart euler = q.ToEulerZYX();
// euler[0] = yaw ≈ 0.5
// euler[1] = pitch ≈ 0.2
// euler[2] = roll ≈ 0.1
```

### Converting to Rotation Matrix
```cpp
auto q = Quaternion::FromAxisAngle(Vec3Cart(0, 1, 0), Constants::PI / 4);

// Get 3×3 rotation matrix
MatrixNM<Real, 3, 3> R = q.ToRotationMatrix();
// R is orthogonal: R^T * R = I, det(R) = 1
```

## Interpolation

### SLERP (Spherical Linear Interpolation)
```cpp
// Two orientations
auto q1 = Quaternion::Identity();
auto q2 = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI);

// Interpolate at t ∈ [0, 1]
auto q_mid = Quaternion::Slerp(q1, q2, 0.5);  // Halfway rotation

// Animation loop
for (Real t = 0; t <= 1.0; t += 0.1) {
    auto q = Quaternion::Slerp(q1, q2, t);
    // Use q for smooth animation...
}
```

SLERP properties:
- **Constant angular velocity** - smooth rotation speed
- **Shortest path** - always takes the short way around
- **Unit quaternion output** - no need to renormalize

### LERP (Linear Interpolation)
```cpp
// Faster but not constant velocity
auto q = Quaternion::Lerp(q1, q2, t);
q.Normalize();  // MUST normalize for rotation use!
```

## Comparison

```cpp
Quaternion p(1, 0, 0, 0);
Quaternion q(1, 0, 0, 0);

bool equal = (p == q);                  // Exact comparison
bool approx = p.IsApprox(q, 1e-10);     // Within tolerance
```

## Output

```cpp
Quaternion q(0.707, 0, 0.707, 0);

std::cout << q << std::endl;
// Output: [0.707, 0i, 0.707j, 0k]

q.Print();  // Same output to std::cout
```

## Practical Examples

### Camera Orientation
```cpp
// Initialize camera looking down -Z axis
Quaternion cameraOrientation = Quaternion::Identity();

// Rotate camera: yaw (look left/right), pitch (look up/down)
void UpdateCamera(Real yawDelta, Real pitchDelta) {
    auto yawRot = Quaternion::FromAxisAngle(Vec3Cart(0, 1, 0), yawDelta);
    auto pitchRot = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), pitchDelta);
    
    // Apply rotations (order matters!)
    cameraOrientation = yawRot * cameraOrientation * pitchRot;
    cameraOrientation.Normalize();
}

// Get camera's forward direction
Vec3Cart GetForward() {
    return cameraOrientation.Rotate(Vec3Cart(0, 0, -1));
}
```

### Smooth Rotation Animation
```cpp
// Animate from start to end orientation over duration
void AnimateRotation(Quaternion start, Quaternion end, Real duration) {
    for (Real t = 0; t <= duration; t += 0.016) {  // ~60fps
        Real s = t / duration;  // Normalized time [0,1]
        
        Quaternion current = Quaternion::Slerp(start, end, s);
        Vec3Cart euler = current.ToEulerZYX();
        
        std::cout << "t=" << t << " yaw=" << euler[0] 
                  << " pitch=" << euler[1] 
                  << " roll=" << euler[2] << "\n";
    }
}
```

### Combining Rotations
```cpp
// Roll, then pitch, then yaw
auto roll = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), 0.1);
auto pitch = Quaternion::FromAxisAngle(Vec3Cart(0, 1, 0), 0.2);
auto yaw = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), 0.3);

// Combined rotation (right-to-left application)
auto combined = yaw * pitch * roll;

// Apply to a vector
Vec3Cart v(1, 0, 0);
Vec3Cart rotated = combined.Rotate(v);
```

## Conventions Summary

| Aspect | Convention |
|--------|------------|
| Storage | `[w, x, y, z]` (scalar first) |
| Multiplication | Hamilton (right-hand rule) |
| Rotation direction | Right-handed (counter-clockwise looking along axis) |
| Angle units | Radians |
| Euler ZYX | Yaw (Z) → Pitch (Y) → Roll (X), intrinsic |

## See Also

- [Vector](Vector.md) - Vec3Cart used for rotation axes
- [MatrixNM](../core/MatrixNM.md) - 3×3 rotation matrix conversion
- [CoordTransf3D](../algorithms/Coordinate_transformations.md) - Full coordinate transformations

## Runnable Examples

See `src/docs_demos/base/docs_demo_quaternions.cpp` for complete working examples.
