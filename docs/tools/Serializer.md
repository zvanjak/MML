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
✅ **Threshold Filtering** - Filter out large vectors in field visualizations  
✅ **Coordinate Systems** - Cartesian coordinates

---

## Highlights
- Unified persistence for functions, fields, curves, particles, and ODE solutions.
- Multiple output formats optimized for plotting and visualization pipelines.
- Detailed `SerializeResult` for robust error handling and diagnostics.
- Precision, thresholds, and sampling controls to tailor output size and quality.
- Interoperable with Tools: Visualizers, ConsolePrinter, and ODE solvers.

---

## Quick Reference Table

| Format | Dimension | Type | Use Case | Method |
|--------|-----------|------|----------|--------|
| REAL_FUNCTION | 1D | Scalar | Single curve plot | `SaveRealFunc()` |
| REAL_FUNCTION_EQUALLY_SPACED | 1D | Scalar | Optimized storage | `SaveRealFuncEquallySpaced()` |
| MULTI_REAL_FUNCTION | 1D | Multiple Scalars | Compare functions | `SaveRealMultiFunc()` |
| PARAMETRIC_CURVE_CARTESIAN_2D | 2D | Curve | Path in 2D | `SaveParamCurveCartesian2D()` |
| PARAMETRIC_CURVE_CARTESIAN_3D | 3D | Curve | Path in 3D | `SaveParamCurveCartesian3D()` |
| SCALAR_FUNCTION_CARTESIAN_2D | 2D | Surface | Height field | `SaveScalarFunc2DCartesian()` |
| VECTOR_FIELD_2D_CARTESIAN | 2D | Vectors | Flow in 2D | `SaveVectorFunc2DCartesian()` |
| VECTOR_FIELD_3D_CARTESIAN | 3D | Vectors | Flow in 3D | `SaveVectorFunc3DCartesian()` |
| PARTICLE_SIMULATION_DATA_2D | 2D | Animation | 2D particles | `SaveParticleSimulation2D()` |
| PARTICLE_SIMULATION_DATA_3D | 3D | Animation | 3D particles | `SaveParticleSimulation3D()` |

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

## Error Handling

### SerializeError Enum

```cpp
enum class SerializeError {
    OK,                   // Operation succeeded
    FILE_NOT_OPENED,      // Failed to open/create file
    INVALID_PARAMETERS,   // Invalid input parameters
    WRITE_FAILED          // Write operation failed
};
```

### SerializeResult Structure

```cpp
struct SerializeResult {
    bool success;             // true if successful
    SerializeError error;     // Error code (if any)
    std::string message;      // Human-readable error message
};
```

### Error Handling Pattern

```cpp
auto result = Serializer::SaveRealFunc(f, "data", 0, 1, 100, "out.txt");
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

### Robust Error Handling

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

## API Reference

### Real Function Serialization

#### SaveRealFunc (Equally-Spaced)

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

#### SaveRealFunc (Arbitrary Points)

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

#### SaveRealFuncEquallySpaced (Compact)

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

**Advantages:**
- ~50% smaller file size (no x-coordinates stored)
- X-coordinates reconstructed from x1, x2, numPoints
- Perfect for visualization tools that support regular grids

**Use Case:** When x-coordinates can be reconstructed from x1, x2, numPoints.

---

### Multi-Function Serialization

Save multiple functions to one file for comparison plots.

#### SaveRealMultiFunc (Function Pointers)

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

#### SaveRealMultiFunc (Interpolated Functions)

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

### Parametric Curves

#### SaveAsParamCurve2D

Save a 2D parametric curve from coordinate vectors.

```cpp
static SerializeResult SaveAsParamCurve2D(
    const Vector<Real>& vec_x,     // X-coordinates
    const Vector<Real>& vec_y,     // Y-coordinates
    std::string title,             // Curve title
    std::string fileName,          // Output file path
    Real t1 = 0.0,                 // Parameter start value
    Real t2 = 1.0                  // Parameter end value
);
```

**File Format:**
```
PARAMETRIC_CURVE_CARTESIAN_2D
<title>
t1: <t1>
t2: <t2>
NumPoints: <n>
<t[0]> <x[0]> <y[0]>
<t[1]> <x[1]> <y[1]>
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

auto result = Serializer::SaveAsParamCurve2D(x, y, "Circle", "circle.txt", 0.0, 2*M_PI);
```

#### SaveParamCurveCartesian2D

Save 2D parametric curve from function.

```cpp
static bool SaveParamCurveCartesian2D(
    const IRealToVectorFunction<2>& f,
    std::string title,
    Real t1, Real t2,
    int numPoints,
    std::string fileName
);
```

**Example (Spiral):**
```cpp
auto spiral = [](Real t) {
    return VectorN<Real, 2>{t * std::cos(t), t * std::sin(t)};
};
RealToVectorFunction<2> curve(spiral);

Serializer::SaveParamCurveCartesian2D(
    curve, "Archimedean Spiral", 0, 4*M_PI, 500, "spiral.txt"
);
```

#### SaveParamCurveCartesian3D

Save 3D parametric curve.

```cpp
static bool SaveParamCurveCartesian3D(
    const IRealToVectorFunction<3>& f,
    std::string title,
    Real t1, Real t2,
    int numPoints,
    std::string fileName
);
```

---

### Scalar Field Serialization

#### SaveScalarFunc2DCartesian

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

**Grid Structure:**
- X varies in outer loop (slower)
- Y varies in inner loop (faster)
- Total points: NumPointsX × NumPointsY

**Example (2D Gaussian):**
```cpp
auto gaussian = [](const Vector<Real>& v) {
    return std::exp(-(v[0]*v[0] + v[1]*v[1]));
};
ScalarFunction<2> f(gaussian);

auto result = Serializer::SaveScalarFunc2DCartesian(
    f, "2D Gaussian", -3, 3, 50, -3, 3, 50, "gaussian2d.txt");
```

#### SaveScalarFunc3DCartesian

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

**Loop Nesting:**
- X: Outer loop (slowest)
- Y: Middle loop
- Z: Inner loop (fastest)
- Total Points: NumPointsX × NumPointsY × NumPointsZ

---

### Vector Field Serialization

#### SaveVectorFunc2DCartesian

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

**With Magnitude Threshold:**

```cpp
static SerializeResult SaveVectorFunc2DCartesian(
    const IVectorFunction<2>& f,
    std::string title,
    Real x1, Real x2, int numPointsX,
    Real y1, Real y2, int numPointsY,
    std::string fileName,
    Real upper_threshold  // Filter out vectors with ||v|| >= threshold
);
```

**File Format:**
```
VECTOR_FIELD_2D_CARTESIAN
<title>
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

**Example with Threshold:**
```cpp
auto electricField = [](const VectorN<Real, 2>& r) {
    Real rMag = r.NormL2();
    if (rMag < 1e-10) return VectorN<Real, 2>{0, 0};
    return r / (rMag * rMag * rMag);  // 1/r² field
};

VectorFunction<2> field(electricField);

// Filter out very large vectors
Serializer::SaveVectorFunc2DCartesian(
    field, "E-field", -5, 5, 30, -5, 5, 30,
    "efield.txt", 50.0  // Exclude ||E|| >= 50
);
```

#### SaveVectorFunc3DCartesian

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

// With threshold variant
static SerializeResult SaveVectorFunc3DCartesian(
    const IVectorFunction<3>& f,
    std::string title,
    Real x1, Real x2, int numPointsX,
    Real y1, Real y2, int numPointsY,
    Real z1, Real z2, int numPointsZ,
    std::string fileName,
    Real upper_threshold
);
```

---

### ODE Solution Serialization

#### SaveODESolutionAsMultiFunc

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

#### SaveODESolutionComponentAsFunc

Save a single component of the solution.

```cpp
static SerializeResult SaveODESolutionComponentAsFunc(
    const ODESystemSolution& sol,
    int compInd,                      // Component index
    std::string compName,             // Component name
    std::string fileName
);
```

**Example:**
```cpp
// After solving ODE system
auto result = Serializer::SaveODESolutionComponentAsFunc(
    sol, 0, "Position vs Time", "position.txt"
);
```

#### SaveODESolAsParametricCurve2D

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

#### SaveODESolAsParametricCurve3D

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

**Example (Lorenz Attractor):**
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

### Particle Simulation Serialization

#### SaveParticleSimulation2D

Save 2D particle animation data.

```cpp
static SerializeResult SaveParticleSimulation2D(
    std::string fileName,
    int numBalls,
    Real width, Real height,
    const std::vector<std::vector<Vector2Cartesian>>& positions,
    const std::vector<std::string>& ballColors,
    const std::vector<Real>& ballRadius,
    Real dT,
    int saveEveryNSteps = 1
);
```

**File Format:**
```
PARTICLE_SIMULATION_DATA_2D
NumBalls: <numBalls>
Width: <width>
Height: <height>
TimeStep: <dT>
NumSteps: <frames>
Step 0 0.000000
<ballIndex> <x> <y>
<ballIndex> <x> <y>
...
Step 1 <dT>
<ballIndex> <x> <y>
<ballIndex> <x> <y>
...
```

**Parameters:**
- `positions[i]`: Position history of particle i
- `ballColors[i]`: Color string (e.g., "red", "#FF0000")
- `ballRadius[i]`: Radius for rendering
- `dT`: Time step between simulation updates
- `saveEveryNSteps`: Frame decimation (default 1 = all frames)

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
    "animation.txt", numBalls, 10.0, 10.0, trajectory,
    ballColors, ballRadii, dT);
```

#### SaveParticleSimulation3D

Save 3D particle animation data.

```cpp
static SerializeResult SaveParticleSimulation3D(
    std::string fileName,
    int numBalls,
    Real width, Real height, Real depth,
    const std::vector<std::vector<Vector3Cartesian>>& positions,
    const std::vector<std::string>& ballColors,
    const std::vector<Real>& ballRadius,
    Real dT,
    int saveEveryNSteps = 1
);
```

---

## Common Patterns

### Comparing Interpolation Methods

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

### Numerical vs Analytical Solution

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

### Vector Field Visualization

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

---

## Usage Examples

### Example 1: Save and Visualize Real Function

```cpp
// Save data
RealFunction sine([](Real x) { return std::sin(x); });
auto result = Serializer::SaveRealFunc(
    sine, "Sine Wave", 0.0, 2*M_PI, 200, "results/sine.txt", 8
);

if (result.success) {
    std::cout << "Data saved successfully!" << std::endl;
}
```

### Example 2: ODE Solution Export (Multiple Views)

```cpp
ODESystemFixedStepSolver<4> solver(system);
auto sol = solver.integrate_RK4(initState, 0.0, 10.0, 0.01);

// Export multiple views
std::vector<std::string> labels = {"x", "y", "vx", "vy"};
Serializer::SaveODESolutionAsMultiFunc(sol, "All", labels, "all.txt");
Serializer::SaveODESolAsParametricCurve2D(sol, "phase.txt", 0, 1, "Phase");
Serializer::SaveODESolutionComponentAsFunc(sol, 0, "X(t)", "x_t.txt");
```

### Example 3: Batch Processing

```cpp
std::vector<IRealFunction*> solutions;
for (int i = 0; i < numCases; ++i) {
    auto sol = solveProblem(params[i]);
    solutions.push_back(&sol);
    
    std::string filename = "case_" + std::to_string(i) + ".txt";
    Serializer::SaveRealFunc(*solutions[i], "Case " + std::to_string(i),
                             0, 10, 100, filename);
}

// Also save comparison
Serializer::SaveRealMultiFunc(solutions, "All Cases", legends,
                              0, 10, 100, "comparison.txt");
```

---

## Integration Patterns

### With Visualizer

The Visualizer class provides static methods that internally use Serializer to save data before launching external visualization tools:

```cpp
// Direct visualization (Visualizer handles serialization internally)
Visualizer::VisualizeRealFunction(f, "sin(x)", 0, 2*M_PI, 100, "sine.txt");

// Or with precision parameter
Visualizer::VisualizeRealFunction(f, "sin(x)", 0, 2*M_PI, 100, "sine.txt", 8);
```

**When to use Serializer directly:**
- Data-only export (no visualization needed)
- Batch processing multiple datasets
- Integration with external tools (Gnuplot, Python, Matlab)
- When visualization will happen in a separate session

**When to use Visualizer:**
- Interactive exploration of results
- Quick visual verification of computations

### With ODE Solvers

```cpp
// Solve system
ODESystemSolution solution = solver.solve(system, init_cond);

// Export in multiple formats for different analysis
std::vector<std::string> legend = {"x", "y", "z"};

// Time series of all components
Serializer::SaveODESolutionAsMultiFunc(
    solution, "State Variables", legend, "states.txt");

// Phase space trajectory
Serializer::SaveODESolAsParametricCurve3D(
    solution, "phase_space.txt", 0, 1, 2, "Phase Portrait");

// Individual component analysis
Serializer::SaveODESolutionComponentAsFunc(
    solution, 0, "X-component", "x_component.txt");
```

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

**Precision Impact:**
- **15 digits**: ~30 bytes/point, scientific accuracy
- **6 digits**: ~18 bytes/point, visualization quality
- **3 digits**: ~15 bytes/point, qualitative plots

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
#include <filesystem>
std::filesystem::create_directories("results/functions");
std::filesystem::create_directories("results/ode");
std::filesystem::create_directories("results/fields");

Serializer::SaveRealFunc(f, title, x1, x2, n, "results/functions/sine.txt");
Serializer::SaveODESolutionAsMultiFunc(sol, title, legend, "results/ode/pendulum.txt");
Serializer::SaveVectorFunc2DCartesian(field, title, x1,x2,nx, y1,y2,ny, "results/fields/gradient.txt");
```

### 6. Use Space-Optimized Format When Appropriate

```cpp
// For regular grids, use equally-spaced format
Serializer::SaveRealFuncEquallySpaced(f, title, x1, x2, n, "data.txt");
// Saves ~33% disk space compared to REAL_FUNCTION
```

### 7. Threshold Filtering for Vector Fields

```cpp
// Avoid visualization clutter near singularities
auto electricField = [](const VectorN<Real, 2>& r) {
    Real rMag = r.NormL2();
    if (rMag < 1e-10) return VectorN<Real, 2>{0, 0};
    return r / (rMag * rMag * rMag);  // 1/r² field
};

VectorFunction<2> field(electricField);

// Filter out very large vectors
Serializer::SaveVectorFunc2DCartesian(
    field, "E-field", -5, 5, 30, -5, 5, 30,
    "efield.txt", 50.0  // Exclude ||E|| >= 50
);
```

### 8. Sampling Density

```cpp
// Match sampling to feature scale
Serializer::SaveRealFunc(smoothFunc, "Smooth", 0, 10, 100, "smooth.txt");     // Sparse OK
Serializer::SaveRealFunc(oscillatory, "Osc", 0, 10, 1000, "oscillatory.txt"); // Dense needed
```

**Sampling Recommendations:**
- **Plotting/visualization**: 100-1000 points usually sufficient
- **Numerical analysis**: Higher resolution (1000-10000 points)
- **2D/3D fields**: Balance resolution vs file size
  - 50×50 = 2500 points ≈ 150KB (reasonable)
  - 200×200 = 40000 points ≈ 2.4MB (getting large)

### 9. Legend Size Matching

```cpp
if (legend.size() != funcs.size()) {
    std::cerr << "Legend count mismatch\n";
    return;
}
auto result = Serializer::SaveRealMultiFunc(funcs, "...", legend, ...);
```

### 10. Frame Decimation for Animations

```cpp
// Save every 5th frame to reduce file size
Serializer::SaveParticleSimulation2D(
    "sim.txt", numBalls, width, height,
    ballPositions, ballColors, ballRadii, dT,
    5  // saveEveryNSteps
);
```

---

## Performance Considerations

### File Size Estimation

| Data Type | Bytes/Point (precision=15) | Bytes/Point (precision=6) |
|-----------|---------------------------|---------------------------|
| Real function | ~40 | ~24 |
| Multi-function (n funcs) | ~(20 + 20*n) | ~(12 + 12*n) |
| 2D scalar field | ~60 | ~36 |
| Vector field 2D | ~80 | ~48 |
| Vector field 3D | ~120 | ~72 |

**Example File Sizes:**
- 1000 points, precision=15: ~30 KB
- 100×100 grid (2D field), precision=15: ~900 KB  
- ODE solution, 10000 steps, 4 vars, precision=15: ~1.2 MB

### Memory Usage

Serializer methods don't load entire datasets into memory - they write incrementally. Safe for very large datasets.

### Write Speed

- I/O bound, typically 1-10 MB/s
- Particle simulations use buffering for efficiency (`ostringstream`)
- For very large datasets, consider reducing precision or sampling

---

## Common Issues

### Issue 1: File Not Created

**Cause:** Directory doesn't exist or no write permissions.

**Solution:**
```cpp
#include <filesystem>
std::filesystem::create_directories("results/subfolder");
auto result = Serializer::SaveRealFunc(f, "...", 0, 1, 100, "results/subfolder/data.txt");
```

### Issue 2: Too Many Data Points

**Cause:** Dense sampling creates huge files.

**Solution:**
```cpp
// Reduce sampling
Serializer::SaveScalarFunc2DCartesian(f, "...", -10, 10, 50, -10, 10, 50, "data.txt");  
// Was 200x200 (40K points), now 50x50 (2500 points)

// Or reduce precision
Serializer::SaveRealFunc(f, "...", 0, 100, 10000, "data.txt", 4);  // 4 digits
```

### Issue 3: Legend Size Mismatch

**Cause:** `legend.size() != funcs.size()`.

**Solution:**
```cpp
if (legend.size() != funcs.size()) {
    std::cerr << "Legend count must match function count\n";
    return;
}
auto result = Serializer::SaveRealMultiFunc(funcs, "...", legend, ...);
```

### Issue 4: Vector Field Singularities

**Cause:** Very large vectors near singularities distort visualization.

**Solution:** Use threshold filtering:
```cpp
Serializer::SaveVectorFunc2DCartesian(
    field, "...", -5, 5, 30, -5, 5, 30,
    "field.txt", 100.0  // Filter ||v|| >= 100
);
```

### Issue 5: Parameter Range Invalid

**Cause:** x1 >= x2 or numPoints < 2.

**Solution:** Validate before calling:
```cpp
if (x1 >= x2 || numPoints < 2) {
    std::cerr << "Invalid parameters\n";
    return;
}
```

---

## File Format Reference

### Format Identifiers

| Identifier | Description |
|-----------|-------------|
| `REAL_FUNCTION` | Single function with x and f(x) |
| `REAL_FUNCTION_EQUALLY_SPACED` | Single function, f(x) only |
| `MULTI_REAL_FUNCTION` | Multiple functions on same domain |
| `PARAMETRIC_CURVE_CARTESIAN_2D` | 2D parametric curve |
| `PARAMETRIC_CURVE_CARTESIAN_3D` | 3D parametric curve |
| `SCALAR_FUNCTION_CARTESIAN_2D` | 2D scalar field |
| `SCALAR_FUNCTION_CARTESIAN_3D` | 3D scalar field |
| `VECTOR_FIELD_2D_CARTESIAN` | 2D vector field (Cartesian) |
| `VECTOR_FIELD_3D_CARTESIAN` | 3D vector field (Cartesian) |
| `PARTICLE_SIMULATION_DATA_2D` | 2D particle animation |
| `PARTICLE_SIMULATION_DATA_3D` | 3D particle animation |

### Common Header Format

Most file formats start with:
```
<TYPE_IDENTIFIER>
<title>
<parameter_info>
<data_rows>
```

All formats are self-describing with metadata in headers.

---

## See Also

- **[SerializerFormats](SerializerFormats.md)** - Detailed file format specifications with examples
- **[Visualizer](Visualizer.md)** - Load and display serialized data
- **[ConsolePrinter](ConsolePrinter.md)** - Format tabular data for console output
- **[ODE Solvers](../algorithms/Differential_equations_solvers.md)** - Generate solutions to serialize
- **[Interpolated Functions](../core/Interpolated_functions.md)** - Functions that can be serialized

---

## Summary

The **Serializer** class provides:

✅ **Comprehensive Coverage** - 23 methods for 12+ data types  
✅ **Robust Error Handling** - SerializeResult pattern with detailed messages  
✅ **Flexible Precision** - Configurable decimal places (1-15+ digits)  
✅ **Cartesian Coordinates** - Consistent coordinate system throughout  
✅ **Threshold Filtering** - Avoid singularities in vector field visualization  
✅ **Space Optimization** - Compact formats for regular grids  
✅ **ODE Integration** - Multiple export formats for phase space analysis  
✅ **Particle Animations** - Frame decimation and metadata support  
✅ **Tool Integration** - Plain text format for Gnuplot, Python, Matlab  
✅ **Incremental Writing** - Memory-efficient for large datasets  

**Key Takeaway:** Professional data persistence for all MML mathematical objects with comprehensive error handling, flexible formatting options, and seamless integration with visualization and analysis tools.
