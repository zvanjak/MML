# Serializer - Data Export and Persistence

## Overview

Static class providing comprehensive data serialization to text files for all MML mathematical objects. Features robust error handling with `SerializeResult` pattern and supports functions, fields, curves, surfaces, ODE solutions, and particle simulations.

**Key Features:**
- Structured error handling (SerializeResult pattern)
- High-precision output (configurable, default 15 digits)
- Multiple data types: real functions, scalar/vector fields, parametric curves, ODE solutions, particles
- Coordinate system support (Cartesian, Spherical)
- Threshold filtering for vector fields

---

## Error Handling

### SerializeResult Structure

```cpp
struct SerializeResult {
    bool success;                 // true if successful
    SerializeError error;         // Error code
    std::string message;          // Human-readable message
};

enum class SerializeError {
    OK,                    // Success
    FILE_NOT_OPENED,       // Cannot open file
    INVALID_PARAMETERS,    // Invalid input
    WRITE_FAILED           // Write operation failed
};
```

**Example:**
```cpp
auto result = Serializer::SaveRealFunc(f, "sin(x)", 0, 2*M_PI, 100, "sine.txt");
if (!result.success) {
    std::cerr << "Error: " << result.message << "\n";
    // Handle error based on result.error code
}
```

---

## Real Functions

### Single Function - Equally Spaced

```cpp
static SerializeResult SaveRealFunc(
    const IRealFunction& f,
    std::string title,
    Real x1, Real x2,
    int numPoints,
    std::string fileName,
    int precision = 15
);
```

**Example:**
```cpp
RealFunction sine([](Real x) { return std::sin(x); });
auto result = Serializer::SaveRealFunc(
    sine, "Sine Wave", 0.0, 2*M_PI, 200, "results/sine.txt", 8
);
```

### Single Function - Custom Points

```cpp
static SerializeResult SaveRealFunc(
    const IRealFunction& f,
    std::string title,
    Vector<Real> points,
    std::string fileName,
    int precision = 15
);
```

### Space-Optimized (Values Only)

```cpp
static SerializeResult SaveRealFuncEquallySpaced(
    const IRealFunction& f,
    std::string title,
    Real x1, Real x2,
    int numPoints,
    std::string fileName,
    int precision = 15
);
```

**Note:** Saves only function values (not x-coordinates). Smaller file size; x-coordinates reconstructed from header.

---

## Multiple Real Functions

Save multiple functions to single file for comparison.

```cpp
static SerializeResult SaveRealMultiFunc(
    const std::vector<IRealFunction*>& funcs,
    std::string title,
    std::vector<std::string> legend,
    Real x1, Real x2,
    int numPoints,
    std::string fileName,
    int precision = 15
);
```

**Specialized overloads:**
- `std::vector<LinearInterpRealFunc>`
- `std::vector<PolynomInterpRealFunc>`
- `std::vector<SplineInterpRealFunc>`

**Example:**
```cpp
std::vector<IRealFunction*> funcs = {&sine, &cosine, &tangent};
std::vector<std::string> legend = {"sin(x)", "cos(x)", "tan(x)"};

auto result = Serializer::SaveRealMultiFunc(
    funcs, "Trigonometric Functions", legend,
    0.0, 2*M_PI, 200, "results/trig.txt"
);
```

---

## Scalar Functions (Surfaces)

### 2D Scalar Field

```cpp
static SerializeResult SaveScalarFunc2DCartesian(
    const IScalarFunction<2>& f,
    std::string title,
    Real x1, Real x2, int numPointsX,
    Real y1, Real y2, int numPointsY,
    std::string fileName
);
```

**Example:**
```cpp
auto gaussian = [](const VectorN<Real, 2>& v) {
    Real x = v[0], y = v[1];
    return std::exp(-(x*x + y*y));
};
ScalarFunction<2> f(gaussian);

auto result = Serializer::SaveScalarFunc2DCartesian(
    f, "2D Gaussian", -3, 3, 50, -3, 3, 50, "gaussian_2d.txt"
);
```

### 3D Scalar Field

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

## Vector Fields

### 2D Vector Field

```cpp
static SerializeResult SaveVectorFunc2DCartesian(
    const IVectorFunction<2>& f,
    std::string title,
    Real x1, Real x2, int numPointsX,
    Real y1, Real y2, int numPointsY,
    std::string fileName
);

// With magnitude threshold
static SerializeResult SaveVectorFunc2DCartesian(
    const IVectorFunction<2>& f,
    std::string title,
    Real x1, Real x2, int numPointsX,
    Real y1, Real y2, int numPointsY,
    std::string fileName,
    Real upper_threshold  // Filter out vectors with ||v|| >= threshold
);
```

**Example:**
```cpp
auto gradientField = [](const VectorN<Real, 2>& v) {
    return VectorN<Real, 2>{v[0], v[1]};  // F(x,y) = (x, y)
};
VectorFunction<2> field(gradientField);

// With threshold to avoid singularities
auto result = Serializer::SaveVectorFunc2DCartesian(
    field, "Gradient Field", -5, 5, 25, -5, 5, 25,
    "field.txt", 10.0  // Exclude vectors with magnitude >= 10
);
```

### 3D Vector Field

```cpp
static SerializeResult SaveVectorFunc3DCartesian(
    const IVectorFunction<3>& f,
    std::string title,
    Real x1, Real x2, int numPointsX,
    Real y1, Real y2, int numPointsY,
    Real z1, Real z2, int numPointsZ,
    std::string fileName
);

// With threshold
static SerializeResult SaveVectorFunc3DCartesian(
    /* ... same parameters ... */
    Real upper_threshold
);
```

### Spherical Vector Field

```cpp
static SerializeResult SaveVectorFuncSpherical(
    const IVectorFunction<3>& f,
    std::string title,
    Real r1, Real r2, int numPointsR,
    Real theta1, Real theta2, int numPointsTheta,
    Real phi1, Real phi2, int numPointsPhi,
    std::string fileName
);

// With threshold
static SerializeResult SaveVectorFuncSpherical(
    /* ... same parameters ... */
    Real upper_threshold
);
```

---

## Parametric Curves

### 2D Parametric Curve

```cpp
static bool SaveParamCurveCartesian2D(
    const IRealToVectorFunction<2>& f,
    std::string title,
    Real t1, Real t2,
    int numPoints,
    std::string fileName
);
```

**Example:**
```cpp
auto spiral = [](Real t) {
    return VectorN<Real, 2>{t * std::cos(t), t * std::sin(t)};
};
RealToVectorFunction<2> curve(spiral);

Serializer::SaveParamCurveCartesian2D(
    curve, "Archimedean Spiral", 0, 4*M_PI, 500, "spiral.txt"
);
```

### 3D Parametric Curve

```cpp
static bool SaveParamCurveCartesian3D(
    const IRealToVectorFunction<3>& f,
    std::string title,
    Real t1, Real t2,
    int numPoints,
    std::string fileName
);
```

### From Component Vectors

```cpp
static SerializeResult SaveAsParamCurve2D(
    const Vector<Real>& vec_x,
    const Vector<Real>& vec_y,
    std::string fileName,
    Real t1 = 0.0,
    Real t2 = 1.0
);
```

**Example:**
```cpp
Vector<Real> xCoords = {0, 1, 2, 3, 4};
Vector<Real> yCoords = {0, 1, 4, 9, 16};

auto result = Serializer::SaveAsParamCurve2D(
    xCoords, yCoords, "parabola.txt", 0.0, 4.0
);
```

---

## ODE System Solutions

### Single Component as Function

```cpp
static SerializeResult SaveODESolutionComponentAsFunc(
    const ODESystemSolution& sol,
    int compInd,  // Component index
    std::string title,
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

### All Components as Multi-Function

```cpp
static SerializeResult SaveODESolutionAsMultiFunc(
    const ODESystemSolution& sol,
    std::string title,
    std::vector<std::string> legend,
    std::string fileName
);
```

**Example:**
```cpp
std::vector<std::string> labels = {"x", "y", "vx", "vy"};
auto result = Serializer::SaveODESolutionAsMultiFunc(
    sol, "State Variables", labels, "states.txt"
);
```

### As 2D Parametric Curve

```cpp
static SerializeResult SaveODESolAsParametricCurve2D(
    const ODESystemSolution& sol,
    std::string fileName,
    int ind1, int ind2,  // Component indices
    std::string title
);
```

**Example:**
```cpp
// Phase portrait: velocity vs position
auto result = Serializer::SaveODESolAsParametricCurve2D(
    sol, "phase_portrait.txt", 0, 1, "Phase Space"
);
```

### As 3D Parametric Curve

```cpp
static SerializeResult SaveODESolAsParametricCurve3D(
    const ODESystemSolution& sol,
    std::string fileName,
    int ind1, int ind2, int ind3,
    std::string title
);
```

---

## Particle Simulations

### 2D Particle Simulation

```cpp
static SerializeResult SaveParticleSimulation2D(
    std::string fileName,
    int numBalls,
    Real width, Real height,
    std::vector<std::vector<Pnt2Cart>> ballPositions,
    std::vector<std::string> ballColors,
    std::vector<Real> ballRadius,
    Real dT,
    int saveEveryNSteps = 1
);
```

### 3D Particle Simulation

```cpp
static SerializeResult SaveParticleSimulation3D(
    std::string fileName,
    int numBalls,
    Real width, Real height, Real depth,
    std::vector<std::vector<Pnt3Cart>> ballPositions,
    std::vector<std::string> ballColors,
    std::vector<Real> ballRadius,
    Real dT,
    int saveEveryNSteps = 1
);
```

---

## Data Format Examples

### Real Function

```
REAL_FUNCTION
Sine Function
x1: 0
x2: 6.28319
NumPoints: 100
0 0
0.0632832 0.0632588
0.126566 0.125987
...
```

### Multi-Function

```
MULTI_REAL_FUNCTION
Comparison
3
sin(x)
cos(x)
tan(x)
x1: 0
x2: 6.28319
NumPoints: 100
0 0 1 0
0.0632832 0.0632588 0.998001 0.0633778
...
```

### Parametric Curve 2D

```
PARAMETRIC_CURVE_CARTESIAN_2D
Circle
t1: 0
t2: 6.28319
NumPoints: 100
0 1 0
0.0632832 0.998001 0.0632588
...
```

### Scalar Field 2D

```
SCALAR_FUNCTION_CARTESIAN_2D
Gaussian
x1: -3
x2: 3
NumPointsX: 50
y1: -3
y2: 3
NumPointsY: 50
-3 -3 0.00012341
-3 -2.87755 0.000152337
...
```

### Vector Field 2D

```
VECTOR_FIELD_2D_CARTESIAN
Electric Field
-5 -5 0.01 0.01
-5 -4.5 0.0123 0.00987
...
(format: x y Fx Fy)
```

---

## Best Practices

### 1. Error Handling

```cpp
auto result = Serializer::SaveRealFunc(f, "data", 0, 1, 100, "out.txt");
if (!result.success) {
    switch (result.error) {
        case Serializer::SerializeError::FILE_NOT_OPENED:
            // Check file permissions, path
            break;
        case Serializer::SerializeError::INVALID_PARAMETERS:
            // Fix parameter ranges
            break;
        case Serializer::SerializeError::WRITE_FAILED:
            // Check disk space
            break;
    }
    std::cerr << result.message << "\n";
}
```

### 2. Precision Selection

```cpp
// High precision for numerical analysis
Serializer::SaveRealFunc(f, "Precise", 0, 1, 100, "data.txt", 15);

// Moderate precision for plotting
Serializer::SaveRealFunc(f, "Plot", 0, 1, 100, "plot.txt", 6);

// Low precision for large datasets
Serializer::SaveRealFunc(f, "Large", 0, 1, 10000, "big.txt", 3);
```

### 3. File Organization

```cpp
// Use subdirectories
Serializer::SaveRealFunc(f, "...", 0, 1, 100, "results/case1/func.txt");
Serializer::SaveScalarFunc2DCartesian(s, "...", -1, 1, 50, -1, 1, 50, "results/case1/field.txt");
```

### 4. Threshold for Vector Fields

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

### 5. Sampling Density

```cpp
// Match sampling to feature scale
Serializer::SaveRealFunc(smoothFunc, "Smooth", 0, 10, 100, "smooth.txt");     // Sparse OK
Serializer::SaveRealFunc(oscillatory, "Osc", 0, 10, 1000, "oscillatory.txt"); // Dense needed
```

---

## Integration Patterns

### With Visualizer

```cpp
// Serialize, then visualize
auto result = Serializer::SaveRealFunc(f, "Function", 0, 10, 200, "data.txt");
if (result.success) {
    Visualizer::VisualizeRealFunction(f, "Function", 0, 10, 200, "data.txt");
}
```

**Note:** Visualizer re-exports data; use Serializer alone for data-only export.

### With ODE Solvers

```cpp
ODESystemFixedStepSolver<4> solver(system);
auto sol = solver.integrate_RK4(initState, 0.0, 10.0, 0.01);

// Export multiple views
Serializer::SaveODESolutionAsMultiFunc(sol, "All", labels, "all.txt");
Serializer::SaveODESolAsParametricCurve2D(sol, "phase.txt", 0, 1, "Phase");
Serializer::SaveODESolutionComponentAsFunc(sol, 0, "X(t)", "x_t.txt");
```

### Batch Processing

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

## Performance Considerations

- **File size**: ~30 bytes/point (precision=15), ~18 bytes/point (precision=6)
- **Write speed**: I/O bound, ~1-10 MB/s typical
- **Memory**: Minimal overhead, streaming writes
- **Large datasets**: Consider reducing precision or sampling

**Example file sizes:**
- 1000 points, precision=15: ~30 KB
- 100x100 grid (2D field), precision=15: ~900 KB  
- ODE solution, 10000 steps, 4 vars, precision=15: ~1.2 MB

---

## Common Issues

### Issue: File Not Created

**Cause:** Directory doesn't exist or no write permissions.

**Solution:**
```cpp
#include <filesystem>
std::filesystem::create_directories("results/subfolder");
auto result = Serializer::SaveRealFunc(f, "...", 0, 1, 100, "results/subfolder/data.txt");
```

### Issue: Too Many Data Points

**Cause:** Dense sampling creates huge files.

**Solution:**
```cpp
// Reduce sampling
Serializer::SaveScalarFunc2DCartesian(f, "...", -10, 10, 50, -10, 10, 50, "data.txt");  
// Was 200x200 (4M points), now 50x50 (2500 points)

// Or reduce precision
Serializer::SaveRealFunc(f, "...", 0, 100, 10000, "data.txt", 4);  // 4 digits
```

### Issue: Legend Size Mismatch

**Cause:** `legend.size() != funcs.size()`.

**Solution:**
```cpp
if (legend.size() != funcs.size()) {
    std::cerr << "Legend count mismatch\n";
    return;
}
auto result = Serializer::SaveRealMultiFunc(funcs, "...", legend, ...);
```

---

## See Also

- **[Visualizer](Visualizer.md)**: Automatically calls Serializer, then launches visualizer
- **[ConsolePrinter](ConsolePrinter.md)**: Formatted terminal output
- **[Timer](Timer_ThreadPool.md)**: Performance measurement

---

## Summary

**Serializer** provides:
- ✅ Robust error handling (SerializeResult pattern)
- ✅ High-precision output (configurable 1-15+ digits)
- ✅ 10+ data types (functions, fields, curves, ODE solutions, particles)
- ✅ Coordinate systems (Cartesian, Spherical)
- ✅ Threshold filtering for vector fields
- ✅ Space optimization (values-only format)
- ✅ Multi-function export (comparison files)
- ✅ ODE solution export (component/multi/phase portrait)
- ✅ Particle simulation export (2D/3D with colors/radii)

**Key takeaway:** Professional data persistence for all MML objects with comprehensive error handling and flexible formatting options.
