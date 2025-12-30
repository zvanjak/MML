# Serializer - Data Persistence for Mathematical Objects

**Header:** `mml/tools/Serializer.h`  
**Namespace:** `MML`  
**Type:** Static utility class

---

## Overview

The Serializer class provides comprehensive functionality for saving mathematical objects (functions, ODE solutions, parametric curves, vector fields) to files in formats suitable for visualization and analysis. All methods are static and return detailed error information through the `SerializeResult` structure.

### Key Features

✅ **Real Functions** - Save 1D functions to file  
✅ **Multi-Functions** - Save multiple functions in one file  
✅ **Scalar Fields** - 2D/3D scalar field data  
✅ **Vector Fields** - 2D/3D vector field visualization data  
✅ **ODE Solutions** - Time-series solutions from differential equations  
✅ **Parametric Curves** - 2D/3D curves from parameters  
✅ **Particle Simulations** - Animation data for particle systems  
✅ **Error Handling** - Detailed success/failure reporting  
✅ **Flexible Precision** - Control decimal places in output  

---

## Quick Start

### Save a Simple Function

```cpp
#include "tools/Serializer.h"
using namespace MML;

// Define a function
auto f = [](Real x) { return std::sin(x); };
RealFunction sine(f);

// Save to file (100 points from 0 to 2π)
auto result = Serializer::SaveRealFunc(sine, "Sine Function", 
                                       0.0, 2*M_PI, 100, 
                                       "sine.txt", 8);

if (!result.success) {
    std::cerr << "Error: " << result.message << std::endl;
}
```

### Save ODE Solution

```cpp
// After solving an ODE system
ODESystemSolution solution = solver.solve(system, initial_conditions);

// Save all components as multi-function
std::vector<std::string> legend = {"x", "y", "vx", "vy"};
auto result = Serializer::SaveODESolutionAsMultiFunc(
    solution, "Pendulum Motion", legend, "pendulum.txt");
```

### Save 2D Parametric Curve

```cpp
Vector<Real> t_vals = {0.0, 0.1, 0.2, 0.3, 0.4};
Vector<Real> x_vals = {0.0, 0.5, 1.0, 1.5, 2.0};
Vector<Real> y_vals = {0.0, 0.25, 1.0, 2.25, 4.0};

auto result = Serializer::SaveAsParamCurve2D(
    x_vals, y_vals, "Parabola", "parabola.txt");
```

---

## API Reference

### Error Handling

#### SerializeError Enum

```cpp
enum class SerializeError {
    OK,                   // Operation succeeded
    FILE_NOT_OPENED,      // Failed to open/create file
    INVALID_PARAMETERS,   // Invalid input parameters
    WRITE_FAILED          // Write operation failed
};
```

#### SerializeResult Structure

```cpp
struct SerializeResult {
    bool success;             // true if successful
    SerializeError error;     // Error code (if any)
    std::string message;      // Human-readable error message
};
```

**Usage Example:**
```cpp
auto result = Serializer::SaveRealFunc(/*...*/);
if (!result.success) {
    switch (result.error) {
        case SerializeError::FILE_NOT_OPENED:
            std::cerr << "Cannot create file: " << result.message << std::endl;
            break;
        case SerializeError::INVALID_PARAMETERS:
            std::cerr << "Bad parameters: " << result.message << std::endl;
            break;
        case SerializeError::WRITE_FAILED:
            std::cerr << "Write error: " << result.message << std::endl;
            break;
    }
}
```

---

## Real Function Serialization

### SaveRealFunc (Equally-Spaced)

Save a real function evaluated at equally-spaced points.

```cpp
static SerializeResult SaveRealFunc(
    const IRealFunction& f,       // Function to serialize
    std::string title,            // Display title for data
    Real x1,                      // Start of range
    Real x2,                      // End of range
    int numPoints,                // Number of evaluation points
    std::string fileName,         // Output file path
    int precision = 15            // Decimal places
);
```

**Preconditions:**
- `x1 < x2`
- `numPoints >= 2`
- `!fileName.empty()`

**File Format:**
```
REAL_FUNCTION
<title>
x1: <x1>
x2: <x2>
NumPoints: <numPoints>
<x0> <f(x0)>
<x1> <f(x1)>
...
```

**Example:**
```cpp
auto sine = [](Real x) { return std::sin(x); };
RealFunction f(sine);

auto result = Serializer::SaveRealFunc(
    f, "sin(x)", 0.0, 2*M_PI, 100, "sine.txt", 8);
```

---

### SaveRealFunc (Arbitrary Points)

Save a real function evaluated at specified (non-equally-spaced) points.

```cpp
static SerializeResult SaveRealFunc(
    const IRealFunction& f,       // Function to serialize
    std::string title,            // Display title
    Vector<Real> points,          // X-coordinates to evaluate
    std::string fileName,         // Output file path
    int precision = 15            // Decimal places
);
```

**Preconditions:**
- `points.size() >= 2`

**Example:**
```cpp
// Evaluate at specific interesting points
Vector<Real> points = {0.0, M_PI/6, M_PI/4, M_PI/3, M_PI/2};
auto result = Serializer::SaveRealFunc(
    f, "Sine at special angles", points, "sine_special.txt", 6);
```

---

### SaveRealFuncEquallySpaced (Compact)

Save only function values (space-optimized, no x-coordinates in output).

```cpp
static SerializeResult SaveRealFuncEquallySpaced(
    const IRealFunction& f,
    std::string title,
    Real x1, Real x2, int numPoints,
    std::string fileName,
    int precision = 15
);
```

**File Format:**
```
REAL_FUNCTION_EQUALLY_SPACED
<title>
x1: <x1>
x2: <x2>
NumPoints: <numPoints>
<f(x0)>
<f(x1)>
...
```

**Use Case:** When x-coordinates can be reconstructed from x1, x2, numPoints.

---

## Multi-Function Serialization

Save multiple functions to one file for comparison plots.

### SaveRealMultiFunc (Function Pointers)

```cpp
static SerializeResult SaveRealMultiFunc(
    const std::vector<IRealFunction*>& funcs,  // Functions to save
    std::string title,                         // Display title
    Real x1, Real x2, int numPoints,           // Range and sampling
    std::vector<std::string> legend,           // Function names
    std::string fileName,
    int precision = 15
);
```

**File Format:**
```
REAL_MULTI_FUNCTION
<title>
NumFunctions: <n>
x1: <x1>
x2: <x2>
NumPoints: <numPoints>
<legend[0]> <legend[1]> ... <legend[n-1]>
<x0> <f0(x0)> <f1(x0)> ... <fn-1(x0)>
<x1> <f0(x1)> <f1(x1)> ... <fn-1(x1)>
...
```

**Example:**
```cpp
auto sine = [](Real x) { return std::sin(x); };
auto cosine = [](Real x) { return std::cos(x); };
RealFunction f1(sine), f2(cosine);

std::vector<IRealFunction*> funcs = {&f1, &f2};
std::vector<std::string> legend = {"sin(x)", "cos(x)"};

auto result = Serializer::SaveRealMultiFunc(
    funcs, "Trig Functions", 0, 2*M_PI, 100, legend, "trig.txt", 8);
```

### SaveRealMultiFunc (Interpolated Functions)

Specialized versions for interpolated function types:

```cpp
// Linear interpolation
static SerializeResult SaveRealMultiFunc(
    const std::vector<LinearInterpRealFunc>& funcs, ...);

// Polynomial interpolation
static SerializeResult SaveRealMultiFunc(
    const std::vector<PolynomInterpRealFunc>& funcs, ...);

// Spline interpolation
static SerializeResult SaveRealMultiFunc(
    const std::vector<SplineInterpRealFunc>& funcs, ...);
```

---

## Parametric Curves

### SaveAsParamCurve2D

Save a 2D parametric curve from coordinate vectors.

```cpp
static SerializeResult SaveAsParamCurve2D(
    const Vector<Real>& vec_x,     // X-coordinates
    const Vector<Real>& vec_y,     // Y-coordinates
    std::string title,
    std::string fileName
);
```

**File Format:**
```
PARAMETRIC_CURVE_2D
<title>
NumPoints: <n>
<x[0]> <y[0]>
<x[1]> <y[1]>
...
```

**Example (Circle):**
```cpp
int n = 100;
Vector<Real> x(n), y(n);
for (int i = 0; i < n; i++) {
    Real t = 2*M_PI * i / n;
    x[i] = std::cos(t);
    y[i] = std::sin(t);
}

auto result = Serializer::SaveAsParamCurve2D(x, y, "Circle", "circle.txt");
```

---

## Scalar Field Serialization

### SaveScalarFunc2DCartesian

Save 2D scalar field data on a Cartesian grid.

```cpp
static SerializeResult SaveScalarFunc2DCartesian(
    const IScalarFunction<2>& f,
    std::string title,
    Real x1, Real x2, int numPointsX,
    Real y1, Real y2, int numPointsY,
    std::string fileName
);
```

**File Format:**
```
SCALAR_FUNC_2D
<title>
x1: <x1> x2: <x2> numX: <numPointsX>
y1: <y1> y2: <y2> numY: <numPointsY>
<x0> <y0> <f(x0,y0)>
<x0> <y1> <f(x0,y1)>
...
```

**Example (2D Gaussian):**
```cpp
auto gaussian = [](const Vector<Real>& v) {
    return std::exp(-(v[0]*v[0] + v[1]*v[1]));
};
ScalarFunction<2> f(gaussian);

auto result = Serializer::SaveScalarFunc2DCartesian(
    f, "2D Gaussian", -3, 3, 50, -3, 3, 50, "gaussian2d.txt");
```

### SaveScalarFunc3DCartesian

Save 3D scalar field data on a Cartesian grid.

```cpp
static SerializeResult SaveScalarFunc3DCartesian(
    const IScalarFunction<3>& f,
    std::string title,
    Real x1, Real x2, int numPointsX,
    Real y1, Real y2, int numPointsY,
    Real z1, Real z2, int numPointsZ,
    std::string fileName
);
```

---

## Vector Field Serialization

### SaveVectorFunc2DCartesian

Save 2D vector field data.

```cpp
static SerializeResult SaveVectorFunc2DCartesian(
    const IVectorFunction<2>& f,
    std::string title,
    Real x1, Real x2, int numPointsX,
    Real y1, Real y2, int numPointsY,
    std::string fileName
);
```

**File Format:**
```
VECTOR_FIELD_2D
<title>
Cartesian
x1: <x1> x2: <x2> numX: <numPointsX>
y1: <y1> y2: <y2> numY: <numPointsY>
<x> <y> <fx(x,y)> <fy(x,y)>
...
```

**Example (Gradient Field):**
```cpp
auto grad_field = [](const Vector<Real>& v) -> Vector<Real> {
    return Vector<Real>{-v[0], -v[1]};  // Points toward origin
};
VectorFunction<2> field(grad_field);

auto result = Serializer::SaveVectorFunc2DCartesian(
    field, "Gradient toward origin", 
    -2, 2, 20, -2, 2, 20, "gradient.txt");
```

### SaveVectorFunc3DCartesian

Save 3D vector field data.

```cpp
static SerializeResult SaveVectorFunc3DCartesian(
    const IVectorFunction<3>& f,
    std::string title,
    Real x1, Real x2, int numPointsX,
    Real y1, Real y2, int numPointsY,
    Real z1, Real z2, int numPointsZ,
    std::string fileName
);
```

---

## ODE Solution Serialization

### SaveODESolutionAsMultiFunc

Save all components of an ODE solution as a multi-function.

```cpp
static SerializeResult SaveODESolutionAsMultiFunc(
    const ODESystemSolution& sol,
    std::string title,
    std::vector<std::string> legend,  // Component names
    std::string fileName
);
```

**Example (Pendulum):**
```cpp
// After solving pendulum ODE
ODESystemSolution solution = solver.solve(pendulum, init_cond);

std::vector<std::string> legend = {"angle", "angular_velocity"};
auto result = Serializer::SaveODESolutionAsMultiFunc(
    solution, "Simple Pendulum", legend, "pendulum.txt");
```

### SaveODESolutionComponentAsFunc

Save a single component of the solution.

```cpp
static SerializeResult SaveODESolutionComponentAsFunc(
    const ODESystemSolution& sol,
    int compInd,                      // Component index
    std::string compName,             // Component name
    std::string fileName
);
```

### SaveODESolAsParametricCurve2D

Save trajectory in phase space (2D).

```cpp
static SerializeResult SaveODESolAsParametricCurve2D(
    const ODESystemSolution& sol,
    std::string fileName,
    int comp_ind_x,                   // Index of x-component
    int comp_ind_y,                   // Index of y-component
    std::string compNameX,
    std::string compNameY
);
```

**Example (Phase Portrait):**
```cpp
// Plot position vs velocity
auto result = Serializer::SaveODESolAsParametricCurve2D(
    solution, "phase_portrait.txt", 
    0, 1,  // Components 0 and 1
    "Position", "Velocity");
```

### SaveODESolAsParametricCurve3D

Save trajectory in 3D space.

```cpp
static SerializeResult SaveODESolAsParametricCurve3D(
    const ODESystemSolution& sol,
    std::string fileName,
    int comp_ind_x, int comp_ind_y, int comp_ind_z,
    std::string compNameX,
    std::string compNameY,
    std::string compNameZ
);
```

---

## Particle Simulation Serialization

### SaveParticleSimulation2D

Save 2D particle animation data.

```cpp
static SerializeResult SaveParticleSimulation2D(
    std::string fileName,
    int numBalls,
    Real width, Real height,
    const std::vector<std::vector<Vector2Cartesian>>& positions
);
```

**File Format:**
```
PARTICLE_SIMULATION_2D
NumParticles: <numBalls>
Width: <width>
Height: <height>
NumFrames: <frames>
# Frame 0
<x0> <y0>
<x1> <y1>
...
# Frame 1
<x0> <y0>
...
```

**Example:**
```cpp
int numBalls = 10;
int numFrames = 100;
std::vector<std::vector<Vector2Cartesian>> trajectory(numFrames);

// Simulate particle motion
for (int frame = 0; frame < numFrames; frame++) {
    trajectory[frame].resize(numBalls);
    for (int i = 0; i < numBalls; i++) {
        // Update particle positions
        trajectory[frame][i] = updatePosition(particles[i], dt);
    }
}

auto result = Serializer::SaveParticleSimulation2D(
    "animation.txt", numBalls, 10.0, 10.0, trajectory);
```

### SaveParticleSimulation3D

Save 3D particle animation data.

```cpp
static SerializeResult SaveParticleSimulation3D(
    std::string fileName,
    int numBalls,
    Real width, Real height, Real depth,
    const std::vector<std::vector<Vector3Cartesian>>& positions
);
```

---

## Usage Examples

### Example 1: Compare Interpolation Methods

```cpp
// Original data
Vector<Real> x_data = {0, 1, 2, 3, 4};
Vector<Real> y_data = {0, 1, 4, 9, 16};

// Create different interpolations
LinearInterpRealFunc linear(x_data, y_data);
PolynomInterpRealFunc poly(x_data, y_data);
SplineInterpRealFunc spline(x_data, y_data);

// Save all for comparison
std::vector<IRealFunction*> funcs = {&linear, &poly, &spline};
std::vector<std::string> legend = {"Linear", "Polynomial", "Spline"};

auto result = Serializer::SaveRealMultiFunc(
    funcs, "Interpolation Comparison", 
    0, 4, 100, legend, "interp_compare.txt", 8);
```

### Example 2: Save Numerical vs Analytical Solution

```cpp
// Analytical solution
auto analytical = [](Real t) { return std::exp(-t); };
RealFunction exact(analytical);

// Numerical solution from ODE solver
ODESystemSolution numerical = solver.solve(ode, init);

// Save both
std::vector<IRealFunction*> funcs = {&exact};
auto result1 = Serializer::SaveRealFunc(
    exact, "Exact", 0, 5, 100, "exact.txt");

auto result2 = Serializer::SaveODESolutionComponentAsFunc(
    numerical, 0, "Numerical", "numerical.txt");
```

### Example 3: Vector Field Visualization

```cpp
// Electric field around two charges
auto electric_field = [](const Vector<Real>& r) -> Vector<Real> {
    Vector<Real> r1{-1, 0}, r2{1, 0};  // Charge positions
    Vector<Real> d1 = r - r1, d2 = r - r2;
    Real dist1 = d1.NormL2(), dist2 = d2.NormL2();
    
    // Field from both charges
    return d1 / (dist1*dist1*dist1) + d2 / (dist2*dist2*dist2);
};
VectorFunction<2> field(electric_field);

auto result = Serializer::SaveVectorFunc2DCartesian(
    field, "Electric Dipole", 
    -3, 3, 30, -3, 3, 30, "dipole.txt");
```

### Example 4: 3D Trajectory

```cpp
// Solve Lorenz system
ODESystemSolution lorenz_sol = solver.solve(lorenz_system, init);

// Save 3D trajectory
auto result = Serializer::SaveODESolAsParametricCurve3D(
    lorenz_sol, "lorenz_attractor.txt",
    0, 1, 2,  // x, y, z components
    "x", "y", "z");
```

---

## File Format Reference

### Common Header Format

Most file formats start with:
```
<TYPE_IDENTIFIER>
<title>
<parameter_info>
<data_rows>
```

### Type Identifiers

- `REAL_FUNCTION` - Single function, x and f(x)
- `REAL_FUNCTION_EQUALLY_SPACED` - Single function, f(x) only
- `REAL_MULTI_FUNCTION` - Multiple functions
- `PARAMETRIC_CURVE_2D` - 2D curve
- `PARAMETRIC_CURVE_3D` - 3D curve
- `SCALAR_FUNC_2D` - 2D scalar field
- `SCALAR_FUNC_3D` - 3D scalar field
- `VECTOR_FIELD_2D` - 2D vector field
- `VECTOR_FIELD_3D` - 3D vector field
- `PARTICLE_SIMULATION_2D` - 2D animation
- `PARTICLE_SIMULATION_3D` - 3D animation

---

## Best Practices

### 1. Always Check Return Value

```cpp
auto result = Serializer::SaveRealFunc(/*...*/);
if (!result.success) {
    std::cerr << "Save failed: " << result.message << std::endl;
    return;
}
std::cout << "Data saved successfully!" << std::endl;
```

### 2. Choose Appropriate Precision

```cpp
// High precision for scientific data
Serializer::SaveRealFunc(f, title, x1, x2, n, file, 15);

// Lower precision for visualization
Serializer::SaveRealFunc(f, title, x1, x2, n, file, 6);

// Very low precision for large datasets
Serializer::SaveRealFunc(f, title, x1, x2, n, file, 3);
```

### 3. Use Descriptive Titles

```cpp
// Good
Serializer::SaveRealFunc(f, "sin(x) from 0 to 2π", ...);

// Better
Serializer::SaveRealFunc(f, "Sine function: amplitude=1, period=2π", ...);
```

### 4. Validate Parameters

```cpp
Real x1 = 0, x2 = 10;
int n = 100;

if (x1 >= x2) {
    std::cerr << "Invalid range: x1 must be < x2" << std::endl;
    return;
}
if (n < 2) {
    std::cerr << "Too few points: need at least 2" << std::endl;
    return;
}

auto result = Serializer::SaveRealFunc(f, title, x1, x2, n, file);
```

### 5. Organize Output Files

```cpp
// Use descriptive filenames and organize in directories
Serializer::SaveRealFunc(f, title, x1, x2, n, "results/functions/sine.txt");
Serializer::SaveODESolutionAsMultiFunc(sol, title, legend, "results/ode/pendulum.txt");
Serializer::SaveVectorFunc2DCartesian(field, title, x1,x2,nx, y1,y2,ny, "results/fields/gradient.txt");
```

---

## Integration with Visualizers

Serialized files are designed to be easily loaded by MML Visualizers:

```cpp
// Save data
Serializer::SaveRealFunc(f, "sin(x)", 0, 2*M_PI, 100, "sine.txt");

// Visualize (in separate application or later)
Visualizer viz;
viz.AddRealFunction("sine.txt");
viz.Show();
```

---

## Error Handling Guide

### Common Errors and Solutions

| Error | Cause | Solution |
|-------|-------|----------|
| `FILE_NOT_OPENED` | Invalid path, permissions, disk full | Check file path, ensure write permissions |
| `INVALID_PARAMETERS` | x1 >= x2, numPoints < 2, empty filename | Validate inputs before calling |
| `WRITE_FAILED` | Disk I/O error, out of space | Check disk space, retry operation |

### Robust Error Handling Pattern

```cpp
auto result = Serializer::SaveRealFunc(f, title, x1, x2, n, fileName);

if (!result.success) {
    std::cerr << "Serialization failed!" << std::endl;
    std::cerr << "Error: " << result.message << std::endl;
    
    // Log error
    logfile << "Failed to save " << fileName << ": " << result.message << std::endl;
    
    // Retry or fallback
    if (result.error == SerializeError::FILE_NOT_OPENED) {
        // Try alternative location
        result = Serializer::SaveRealFunc(f, title, x1, x2, n, "backup/" + fileName);
    }
}
```

---

## Performance Considerations

### File Size Estimation

- **Real function**: ~40 bytes per point (with precision=15)
- **Multi-function (n funcs)**: ~(20 + 20*n) bytes per point
- **2D scalar field (nx×ny)**: ~60 bytes per grid point
- **Vector field 2D (nx×ny)**: ~80 bytes per grid point

### Memory Usage

Serializer methods don't load entire datasets into memory - they write incrementally. Safe for very large datasets.

### Sampling Recommendations

```cpp
// For plotting/visualization: 100-1000 points usually sufficient
Serializer::SaveRealFunc(f, title, x1, x2, 500, file);

// For numerical analysis: higher resolution
Serializer::SaveRealFunc(f, title, x1, x2, 10000, file);

// For 2D/3D fields: balance resolution vs file size
// 50×50 = 2500 points ≈ 150KB (reasonable)
// 200×200 = 40000 points ≈ 2.4MB (getting large)
```

---

## See Also

- **Visualizer** - Load and display serialized data
- **ConsolePrinter** - Format tabular data for console output
- **ODESolver** - Generate solutions to serialize
- **InterpolatedFunction** - Functions that can be serialized

---

## Version History

**Current Version**
- Comprehensive support for 1D, 2D, 3D functions
- ODE solution serialization
- Particle simulation data
- Detailed error reporting with SerializeResult
- Flexible precision control
