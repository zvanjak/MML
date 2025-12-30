# MML Serialization Basics

## Overview

The `Serializer` class provides a comprehensive suite of static methods for exporting mathematical functions, curves, and simulation data to human-readable text files. All formats are designed for easy parsing by visualization tools and data analysis scripts.

## File Format Summary

### 1. **REAL_FUNCTION** - 1D Real Functions
**Use Case:** Saving single real-valued function values

**Format:**
```
REAL_FUNCTION
<title>
x1: <start_x>
x2: <end_x>
NumPoints: <number_of_points>
x1_value y1_value
x2_value y2_value
...
```

**Key Points:**
- Stores both x and y coordinates
- Used for any 1D function: f(x) → y
- Points can be equally-spaced or arbitrary

**Methods:**
- `SaveRealFunc()` - From equally-spaced interval or arbitrary point list

---

### 2. **REAL_FUNCTION_EQUALLY_SPACED** - Optimized 1D Functions
**Use Case:** Efficient storage when points are equally distributed

**Format:**
```
REAL_FUNCTION_EQUALLY_SPACED
<title>
x1: <start_x>
x2: <end_x>
NumPoints: <number_of_points>
y1_value
y2_value
...
```

**Key Points:**
- Only stores y-values; x-values are reconstructed from x1, x2, numPoints
- More compact than REAL_FUNCTION
- Perfect for visualization tools that support regular grids

**Methods:**
- `SaveRealFuncEquallySpaced()`

---

### 3. **MULTI_REAL_FUNCTION** - Multiple 1D Functions
**Use Case:** Comparing multiple functions on same domain

**Format:**
```
MULTI_REAL_FUNCTION
<title>
<num_functions>
<legend_1>
<legend_2>
...
x1: <start_x>
x2: <end_x>
NumPoints: <number_of_points>
x_value y1_value y2_value ... yn_value
x_value y1_value y2_value ... yn_value
...
```

**Key Points:**
- All functions share same x-axis domain
- Includes legend names for each function
- Supports interpolation objects (LinearInterpRealFunc, PolynomInterpRealFunc, SplineInterpRealFunc)

**Methods:**
- `SaveRealMultiFunc()` - with function pointers or interpolation objects

---

### 4. **PARAMETRIC_CURVE_CARTESIAN_2D** - 2D Parametric Curves
**Use Case:** Curves defined as (x(t), y(t))

**Format:**
```
PARAMETRIC_CURVE_CARTESIAN_2D
<title>
t1: <start_t>
t2: <end_t>
NumPoints: <number_of_points>
t1_value x1_value y1_value
t2_value x2_value y2_value
...
```

**Key Points:**
- Parameter t goes from t1 to t2
- Each point is (t, x(t), y(t))
- Ideal for path visualization in 2D

**Methods:**
- `SaveParamCurve<2>()` / `SaveParamCurveCartesian2D()`
- `SaveAsParamCurve()` - from pre-computed point vectors

---

### 5. **PARAMETRIC_CURVE_CARTESIAN_3D** - 3D Parametric Curves
**Use Case:** Curves defined as (x(t), y(t), z(t))

**Format:**
```
PARAMETRIC_CURVE_CARTESIAN_3D
<title>
t1: <start_t>
t2: <end_t>
NumPoints: <number_of_points>
t_value x_value y_value z_value
t_value x_value y_value z_value
...
```

**Key Points:**
- Parameter t goes from t1 to t2
- Each point is (t, x(t), y(t), z(t))
- Ideal for 3D path visualization

**Methods:**
- `SaveParamCurve<3>()` / `SaveParamCurveCartesian3D()`

---

### 6. **SCALAR_FUNCTION_CARTESIAN_2D** - 2D Scalar Fields
**Use Case:** Surface visualization: f(x,y) → z

**Format:**
```
SCALAR_FUNCTION_CARTESIAN_2D
<title>
x1: <start_x>
x2: <end_x>
NumPointsX: <num_x_points>
y1: <start_y>
y2: <end_y>
NumPointsY: <num_y_points>
x1 y1 z_value
x1 y2 z_value
...
x2 y2 z_value
```

**Key Points:**
- Represents surfaces in 2D domain
- Z-axis loop varies faster (inner loop)
- Used for contour plots and 3D surface visualization

**Methods:**
- `SaveScalarFunc2DCartesian()`

---

### 7. **SCALAR_FUNCTION_CARTESIAN_3D** - 3D Scalar Fields
**Use Case:** Volume scalar field: f(x,y,z) → value

**Format:**
```
SCALAR_FUNCTION_CARTESIAN_3D
<title>
x1: <start_x>
x2: <end_x>
NumPointsX: <num_x_points>
y1: <start_y>
y2: <end_y>
NumPointsY: <num_y_points>
z1: <start_z>
z2: <end_z>
NumPointsZ: <num_z_points>
x y z value
x y z value
...
```

**Key Points:**
- 3D volume data representation
- Triple nested loop (z varies fastest)
- Used for volume visualization or iso-surface extraction

**Methods:**
- `SaveScalarFunc3DCartesian()`

---

### 8. **VECTOR_FIELD_2D_CARTESIAN** - 2D Vector Fields
**Use Case:** Arrows/vectors at each point: F(x,y) → (vx, vy)

**Format:**
```
VECTOR_FIELD_2D_CARTESIAN
<title>
x1 y1 vx1 vy1
x1 y2 vx2 vy2
...
x2 y2 vxn vyn
```

**Key Points:**
- Each point has associated 2D vector
- Optional threshold filtering (only include vectors with magnitude < threshold)
- Used for flow field visualization

**Methods:**
- `SaveVectorFunc2D()` / `SaveVectorFunc2DCartesian()`
- Both with and without upper_threshold variants

---

### 9. **VECTOR_FIELD_3D_CARTESIAN** - 3D Vector Fields
**Use Case:** 3D vector fields: F(x,y,z) → (vx, vy, vz)

**Format:**
```
VECTOR_FIELD_3D_CARTESIAN
<title>
x y z vx vy vz
x y z vx vy vz
...
```

**Key Points:**
- Each point has associated 3D vector
- Optional threshold filtering
- Used for 3D flow visualization

**Methods:**
- `SaveVectorFunc3D()` / `SaveVectorFunc3DCartesian()`
- Both with and without upper_threshold variants

---

### 10. **VECTOR_FIELD_SPHERICAL** - Spherical Coordinate Vector Fields
**Use Case:** Vector fields in spherical coordinates: F(r,θ,φ) → (vr, vθ, vφ)

**Format:**
```
VECTOR_FIELD_SPHERICAL
<title>
r θ φ vr vθ vφ
r θ φ vr vθ vφ
...
```

**Key Points:**
- Input coordinates: r (radius), θ (azimuth), φ (polar angle)
- Vector components in spherical basis
- Useful for physics/electromagnetism problems

**Methods:**
- `SaveVectorFuncSpherical()` - with and without threshold

---

### 11. **ODE Solution Formats**

#### Component as REAL_FUNCTION
Save single ODE solution component as if it were a 1D function:
```
REAL_FUNCTION
<title>
x1: <t_start>
x2: <t_end>
NumPoints: <num_steps>
t1 x1(t1)
t2 x1(t2)
...
```

**Methods:** `SaveODESolutionComponentAsFunc()`

#### All Components as MULTI_REAL_FUNCTION
Save entire ODE solution:
```
MULTI_REAL_FUNCTION
<title>
<dimension>
<component_1_name>
<component_2_name>
...
x1: <t_start>
x2: <t_end>
NumPoints: <num_steps>
t x1(t) x2(t) ... xn(t)
t x1(t) x2(t) ... xn(t)
...
```

**Methods:** `SaveODESolutionAsMultiFunc()`

#### Phase Space 2D
Project two components into 2D parametric space:
```
PARAMETRIC_CURVE_CARTESIAN_2D
<title>
t1: <t_start>
t2: <t_end>
NumPoints: <num_steps>
t x_comp1(t) x_comp2(t)
...
```

**Methods:** `SaveODESolAsParametricCurve2D()`

#### Phase Space 3D
Project three components into 3D parametric space:
```
PARAMETRIC_CURVE_CARTESIAN_3D
<title>
t1: <t_start>
t2: <t_end>
NumPoints: <num_steps>
t x_comp1(t) x_comp2(t) x_comp3(t)
...
```

**Methods:** `SaveODESolAsParametricCurve3D()`

---

### 12. **PARTICLE_SIMULATION_DATA_2D** - 2D Particle System
**Use Case:** Animation frames of moving particles

**Format:**
```
PARTICLE_SIMULATION_DATA_2D
Width: <width>
Height: <height>
NumBalls: <count>
Ball_1 <color> <radius>
Ball_2 <color> <radius>
...
NumSteps: <frame_count>
Step 0 0.0
0 x0 y0
1 x1 y1
...
Step 1 0.01
0 x0' y0'
...
```

**Key Points:**
- Includes particle metadata (color, radius)
- Each step records all particle positions
- `saveEveryNSteps` parameter for frame skipping
- Container dimensions define bounds

**Methods:**
- `SaveParticleSimulation2D()`

---

### 13. **PARTICLE_SIMULATION_DATA_3D** - 3D Particle System
**Use Case:** Animation frames in 3D space

**Format:**
```
PARTICLE_SIMULATION_DATA_3D
Width: <width>
Height: <height>
Depth: <depth>
NumBalls: <count>
Ball_1 <color> <radius>
Ball_2 <color> <radius>
...
NumSteps: <frame_count>
Step 0 0.0
0 x0 y0 z0
1 x1 y1 z1
...
Step 1 0.01
...
```

**Key Points:**
- 3D equivalent of 2D particle system
- Includes depth dimension
- Used for collision simulations, molecular dynamics, etc.

**Methods:**
- `SaveParticleSimulation3D()`

---

## Quick Reference Table

| Format | Dimension | Type | Use Case |
|--------|-----------|------|----------|
| REAL_FUNCTION | 1D | Scalar | Single curve |
| REAL_FUNCTION_EQUALLY_SPACED | 1D | Scalar | Optimized storage |
| MULTI_REAL_FUNCTION | 1D | Multiple Scalars | Compare curves |
| PARAMETRIC_CURVE_CARTESIAN_2D | 2D | Curve | Path in 2D |
| PARAMETRIC_CURVE_CARTESIAN_3D | 3D | Curve | Path in 3D |
| SCALAR_FUNCTION_CARTESIAN_2D | 2D | Surface | Height field |
| SCALAR_FUNCTION_CARTESIAN_3D | 3D | Volume | Volumetric data |
| VECTOR_FIELD_2D_CARTESIAN | 2D | Vectors | Flow in 2D |
| VECTOR_FIELD_3D_CARTESIAN | 3D | Vectors | Flow in 3D |
| VECTOR_FIELD_SPHERICAL | Spherical | Vectors | Physics fields |
| PARTICLE_SIMULATION_DATA_2D | 2D | Animation | 2D particles |
| PARTICLE_SIMULATION_DATA_3D | 3D | Animation | 3D particles |

---

## Common Patterns

### Saving a Function
```cpp
// 1D function
Serializer::SaveRealFunc(myFunc, "My Function", -5.0, 5.0, 100, "output.txt");

// Multiple 1D functions
std::vector<IRealFunction*> funcs = {&f1, &f2, &f3};
Serializer::SaveRealMultiFunc(funcs, "Comparison", {"f(x)", "g(x)", "h(x)"}, 
                              -5.0, 5.0, 100, "output.txt");
```

### Saving a 2D Curve
```cpp
Serializer::SaveParamCurveCartesian2D(myCurve, "Spiral", 0.0, 10.0, 200, "spiral.txt");
```

### Saving a 3D Surface
```cpp
Serializer::SaveScalarFunc3DCartesian(myScalar, "Potential Field",
                                     -5, 5, 50,
                                     -5, 5, 50,
                                     -5, 5, 50,
                                     "field.txt");
```

### Saving a Vector Field
```cpp
// 2D field with magnitude filtering
Serializer::SaveVectorFunc2DCartesian(myField, "Gradient",
                                     -3, 3, 30,
                                     -3, 3, 30,
                                     "field.txt", 10.0);  // max magnitude: 10.0
```

### Saving ODE Solution
```cpp
Serializer::SaveODESolutionAsMultiFunc(solution, "ODE System", 
                                      {"x(t)", "y(t)", "z(t)"}, "solution.txt");

// Phase space projection
Serializer::SaveODESolAsParametricCurve3D(solution, "phase_space.txt", 0, 1, 2, "Phase Space");
```

### Saving Particle Simulation
```cpp
Serializer::SaveParticleSimulation2D("sim.txt", numBalls, width, height,
                                    ballPositions, ballColors, ballRadii, dT);
```

---

## Notes

- All methods return `bool` indicating success/failure
- File format headers enable automatic tool detection
- Text format allows easy integration with plotting tools (Gnuplot, Matlab, Python, etc.)
- Floating-point precision follows system default `Real` type precision
