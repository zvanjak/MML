# Visualizer - Data Visualization and Export

## Overview

The `Visualizer` class provides a high-level interface for visualizing mathematical objects and simulation results by exporting data to files and launching external visualization applications. It acts as a bridge between MML's computational components and platform-specific visualization tools.

**Key Features:**
- **Multiple visualization types**: Real functions, scalar/vector fields, parametric curves, surfaces, particles, ODE solutions
- **External visualizer integration**: Launches WPF (Windows) or FLTK/Qt (Linux) applications
- **Data export**: Serializes data to text files for external processing
- **Multi-object visualization**: Plot multiple functions/curves simultaneously
- **Cross-platform path management**: Automatic path resolution for different platforms

**Architecture:**
```
Visualizer (static class)
├── Data Serialization (via Serializer)
├── Path Management (platform-specific)
└── External Process Launch (system calls)
```

---

## Core Concepts

### Visualization Workflow

1. **Serialize**: Convert mathematical object to text file using `Serializer`
2. **Launch**: Execute external visualizer application with file path
3. **Display**: Visualizer app reads data and renders graphics

```cpp
// Example: Visualize a real function
IRealFunction& f = /* ... */;
Visualizer::VisualizeRealFunction(f, "sin(x)", 0.0, 2*M_PI, 100, "sine.txt");
// Creates: results/sine.txt
// Launches: visualizers/RealFuncVisualizer.exe "results/sine.txt"
```

### Path Configuration

The Visualizer uses getter functions for cross-platform path resolution:

```cpp
static std::string _pathResultFiles;           // results/
static std::string _pathRealFuncViz;            // path to RealFuncVisualizer
static std::string _pathSurfaceViz;             // path to SurfaceVisualizer
static std::string _pathParametricCurve3DViz;   // path to 3D curve visualizer
static std::string _pathParametricCurve2DViz;   // path to 2D curve visualizer
static std::string _pathVectorField2DViz;       // path to 2D vector field viz
static std::string _pathVectorField3DViz;       // path to 3D vector field viz
static std::string _pathParticle2DViz;          // path to 2D particle sim viz
static std::string _pathParticle3DViz;          // path to 3D particle sim viz
```

**Note:** Path getters are defined in `MMLBase.h` and detect platform (Windows/Linux) automatically.

---

## API Reference

### Real Functions - Single Function

Visualize functions f: ℝ → ℝ.

#### Equally-Spaced Sampling

```cpp
static void VisualizeRealFunction(
    const IRealFunction& f,
    std::string title,
    Real x1,              // Start of interval
    Real x2,              // End of interval
    int numPoints,        // Number of sample points
    std::string fileName  // Output filename (in results/)
);
```

**Example:**
```cpp
// Visualize sin(x) from 0 to 2π with 100 points
auto sine = [](Real x) { return std::sin(x); };
RealFunction f(sine);

Visualizer::VisualizeRealFunction(
    f, 
    "Sine Function", 
    0.0, 2*M_PI, 
    100, 
    "sine_curve.txt"
);
```

#### Custom Sampling Points

```cpp
static void VisualizeRealFunction(
    const IRealFunction& f,
    std::string title,
    Vector<Real> points,  // Custom x-coordinates
    std::string fileName
);
```

**Example:**
```cpp
// Visualize at specific points
Vector<Real> xValues = {0.0, 0.5, 1.0, 1.5, 2.0, 3.14159};
Visualizer::VisualizeRealFunction(f, "Sample Points", xValues, "samples.txt");
```

---

### Real Functions - Multiple Functions

Visualize multiple functions on the same plot for comparison.

#### Multiple Interface-Based Functions

```cpp
static void VisualizeMultiRealFunction(
    std::vector<IRealFunction*> funcs,
    std::string title,
    std::vector<std::string> func_legend,  // Legend labels
    Real x1, Real x2, 
    int numPoints,
    std::string fileName
);
```

**Example:**
```cpp
// Compare sin, cos, tan
RealFunction sine([](Real x) { return std::sin(x); });
RealFunction cosine([](Real x) { return std::cos(x); });
RealFunction tangent([](Real x) { return std::tan(x); });

std::vector<IRealFunction*> funcs = {&sine, &cosine, &tangent};
std::vector<std::string> legend = {"sin(x)", "cos(x)", "tan(x)"};

Visualizer::VisualizeMultiRealFunction(
    funcs, 
    "Trigonometric Functions", 
    legend, 
    0.0, 2*M_PI, 
    200, 
    "trig_comparison.txt"
);
```

#### Multiple Interpolated Functions

Specialized overloads for interpolated function types:

```cpp
// Linear interpolation functions
static void VisualizeMultiRealFunction(
    std::vector<LinearInterpRealFunc> funcs,
    std::string title,
    std::vector<std::string> func_legend,
    Real x1, Real x2, int numPoints,
    std::string fileName
);

// Polynomial interpolation functions
static void VisualizeMultiRealFunction(
    std::vector<PolynomInterpRealFunc> funcs,
    std::string title,
    std::vector<std::string> func_legend,
    Real x1, Real x2, int numPoints,
    std::string fileName
);

// Spline interpolation functions
static void VisualizeMultiRealFunction(
    std::vector<SplineInterpRealFunc> funcs,
    std::string title,
    std::vector<std::string> func_legend,
    Real x1, Real x2, int numPoints,
    std::string fileName
);
```

**Example:**
```cpp
// Visualize different interpolation methods
Vector<Real> xData = {0, 1, 2, 3, 4};
Vector<Real> yData = {0, 1, 4, 9, 16};

LinearInterpRealFunc linear(xData, yData);
PolynomInterpRealFunc poly(xData, yData);
SplineInterpRealFunc spline(xData, yData);

std::vector<SplineInterpRealFunc> funcs = {spline};
std::vector<std::string> legend = {"Cubic Spline"};

Visualizer::VisualizeMultiRealFunction(
    funcs, 
    "Interpolation Methods", 
    legend, 
    0.0, 4.0, 
    100, 
    "interp_comparison.txt"
);
```

#### Separate File Export (Advanced)

Exports each function to a separate file and visualizes all together:

```cpp
static void VisualizeMultiRealFunctionSeparately(
    std::vector<LinearInterpRealFunc> funcs,
    std::string title,
    std::vector<std::string> func_legend,
    Real x1, Real x2, int numPoints,
    std::string fileName  // Base name; "_1.txt", "_2.txt" appended
);
```

**Use Case:** When functions have different valid domains or require separate processing.

**Example:**
```cpp
// Each function saved to separate file
Visualizer::VisualizeMultiRealFunctionSeparately(
    funcs, 
    "Separate Export", 
    legend, 
    0.0, 10.0, 
    100, 
    "func"
);
// Creates: results/func_1.txt, results/func_2.txt, results/func_3.txt
// Launches visualizer with all three files
```

---

### Scalar Functions (Surfaces)

Visualize scalar fields f: ℝ² → ℝ as 3D surfaces.

```cpp
static void VisualizeScalarFunc2DCartesian(
    const IScalarFunction<2>& func,
    std::string title,
    Real x1, Real x2, int numPointsX,  // X-axis range and sampling
    Real y1, Real y2, int numPointsY,  // Y-axis range and sampling
    std::string fileName
);
```

**Example:**
```cpp
// Visualize f(x,y) = sin(x) * cos(y)
auto surfFunc = [](const VectorN<Real, 2>& v) {
    return std::sin(v[0]) * std::cos(v[1]);
};
ScalarFunction<2> surface(surfFunc);

Visualizer::VisualizeScalarFunc2DCartesian(
    surface,
    "sin(x)*cos(y)",
    -M_PI, M_PI, 50,     // X: [-π, π] with 50 points
    -M_PI, M_PI, 50,     // Y: [-π, π] with 50 points
    "wave_surface.txt"
);
```

**Output Format:** Grid of (x, y, f(x,y)) values for surface plotting.

---

### Vector Fields

Visualize vector fields F: ℝⁿ → ℝⁿ as arrow plots.

#### 2D Vector Fields

```cpp
static void VisualizeVectorField2DCartesian(
    const IVectorFunction<2>& func,
    std::string title,
    Real x1, Real x2, int numPointsX,
    Real y1, Real y2, int numPointsY,
    std::string fileName
);
```

**Example:**
```cpp
// Visualize gradient field: F(x,y) = (x, y)
auto gradField = [](const VectorN<Real, 2>& v) {
    return VectorN<Real, 2>{v[0], v[1]};
};
VectorFunction<2> field(gradField);

Visualizer::VisualizeVectorField2DCartesian(
    field,
    "Radial Field",
    -5.0, 5.0, 20,
    -5.0, 5.0, 20,
    "radial_field.txt"
);
```

**Visualization:** Arrows at each grid point showing direction and magnitude.

#### 3D Vector Fields

```cpp
static void VisualizeVectorField3DCartesian(
    const IVectorFunction<3>& func,
    std::string title,
    Real x1, Real x2, int numPointsX,
    Real y1, Real y2, int numPointsY,
    Real z1, Real z2, int numPointsZ,
    std::string fileName
);
```

**Example:**
```cpp
// Magnetic field around current loop
auto magneticField = [](const VectorN<Real, 3>& r) {
    // Simplified B-field calculation
    Real x = r[0], y = r[1], z = r[2];
    Real rho = std::sqrt(x*x + y*y);
    return VectorN<Real, 3>{-y/rho, x/rho, 0.0};
};
VectorFunction<3> bField(magneticField);

Visualizer::VisualizeVectorField3DCartesian(
    bField,
    "Magnetic Field",
    -2.0, 2.0, 10,
    -2.0, 2.0, 10,
    -1.0, 1.0, 10,
    "magnetic_field.txt"
);
```

---

### Parametric Curves

Visualize curves in 2D and 3D space.

#### 2D Parametric Curves - Single Curve

```cpp
static void VisualizeParamCurve2D(
    const IRealToVectorFunction<2>& f,
    std::string title,
    Real t1, Real t2,      // Parameter range
    int numPoints,
    std::string fileName
);
```

**Example:**
```cpp
// Circle: r(t) = (cos(t), sin(t))
auto circle = [](Real t) {
    return VectorN<Real, 2>{std::cos(t), std::sin(t)};
};
RealToVectorFunction<2> curve(circle);

Visualizer::VisualizeParamCurve2D(
    curve,
    "Unit Circle",
    0.0, 2*M_PI,
    100,
    "circle.txt"
);
```

#### 2D Parametric Curves - Multiple Curves

```cpp
static void VisualizeMultiParamCurve2D(
    std::vector<IRealToVectorFunction<2>*> curves,
    std::string title,
    Real t1, Real t2,
    int numPoints,
    std::string fileName
);

// With custom legend
static void VisualizeMultiParamCurve2D(
    std::vector<IRealToVectorFunction<2>*> curves,
    std::vector<std::string> legend,
    Real t1, Real t2,
    int numPoints,
    std::string fileName
);

// From pre-saved files
static void VisualizeMultiParamCurve2D(
    std::vector<std::string> fileNames  // Relative to results/
);
```

**Example:**
```cpp
// Compare different Lissajous curves
auto lissajous1 = [](Real t) {
    return VectorN<Real, 2>{std::sin(t), std::sin(2*t)};
};
auto lissajous2 = [](Real t) {
    return VectorN<Real, 2>{std::sin(3*t), std::sin(4*t)};
};

RealToVectorFunction<2> curve1(lissajous1);
RealToVectorFunction<2> curve2(lissajous2);

std::vector<IRealToVectorFunction<2>*> curves = {&curve1, &curve2};
std::vector<std::string> legend = {"1:2 Ratio", "3:4 Ratio"};

Visualizer::VisualizeMultiParamCurve2D(
    curves,
    legend,
    0.0, 2*M_PI,
    200,
    "lissajous"
);
```

#### 3D Parametric Curves - Single Curve

```cpp
static void VisualizeParamCurve3D(
    const IRealToVectorFunction<3>& f,
    std::string title,
    Real t1, Real t2,
    int numPoints,
    std::string fileName
);
```

**Example:**
```cpp
// Helix: r(t) = (cos(t), sin(t), t)
auto helix = [](Real t) {
    return VectorN<Real, 3>{std::cos(t), std::sin(t), t};
};
RealToVectorFunction<3> curve(helix);

Visualizer::VisualizeParamCurve3D(
    curve,
    "Helix",
    0.0, 4*M_PI,
    200,
    "helix.txt"
);
```

#### 3D Parametric Curves - Multiple Curves

```cpp
static void VisualizeMultiParamCurve3D(
    std::vector<IRealToVectorFunction<3>*> curves,
    std::string title,
    Real t1, Real t2,
    int numPoints,
    std::string fileName
);

// From pre-saved files
static void VisualizeMultiParamCurve3D(
    std::vector<std::string> fileNames
);
```

**Example:**
```cpp
// Visualize particle trajectories
std::vector<IRealToVectorFunction<3>*> trajectories = {&path1, &path2, &path3};

Visualizer::VisualizeMultiParamCurve3D(
    trajectories,
    "Particle Paths",
    0.0, 10.0,
    500,
    "trajectories"
);
```

---

### ODE System Solutions

Specialized visualizations for ODE system solutions.

#### Single Component as Real Function

```cpp
static void VisualizeODESysSolCompAsFunc(
    const ODESystemSolution& sol,
    int compInd,           // Component index (0-based)
    std::string title,
    std::string fileName
);
```

**Example:**
```cpp
// Solve and visualize position (component 0) vs time
ODESystemSolution sol = /* solve ODE system */;

Visualizer::VisualizeODESysSolCompAsFunc(
    sol,
    0,  // Position component
    "Position vs Time",
    "position_time.txt"
);
```

#### All Components as Multi-Function

```cpp
static void VisualizeODESysSolAsMultiFunc(
    const ODESystemSolution& sol,
    std::string title,
    std::vector<std::string> legend,  // Labels for each component
    std::string fileName
);
```

**Example:**
```cpp
// Visualize all state variables
std::vector<std::string> labels = {"x", "y", "vx", "vy"};

Visualizer::VisualizeODESysSolAsMultiFunc(
    sol,
    "State Variables",
    labels,
    "all_states.txt"
);
```

#### Phase Space - 2D Parametric Curve

```cpp
static void VisualizeODESysSolAsParamCurve2(
    const ODESystemSolution& sol,
    int ind1, int ind2,    // Component indices for x and y
    std::string title,
    std::string fileName
);
```

**Example:**
```cpp
// Phase portrait: velocity vs position
Visualizer::VisualizeODESysSolAsParamCurve2(
    sol,
    0,  // Position (x-axis)
    2,  // Velocity (y-axis)
    "Phase Portrait",
    "phase_portrait.txt"
);
```

**Platform Note:** Currently Windows-only. Linux implementation pending.

#### Phase Space - 3D Parametric Curve

```cpp
static void VisualizeODESysSolAsParamCurve3(
    const ODESystemSolution& sol,
    int ind1, int ind2, int ind3,  // Component indices for x, y, z
    std::string title,
    std::string fileName
);
```

**Example:**
```cpp
// 3D trajectory from ODE solution
Visualizer::VisualizeODESysSolAsParamCurve3(
    sol,
    0, 1, 2,  // x, y, z components
    "3D Trajectory",
    "trajectory_3d.txt"
);
```

**Platform Note:** Currently Windows-only. Linux implementation pending.

---

### Particle Simulations

Visualize pre-computed particle simulation data.

#### 2D Particle Simulation

```cpp
static void VisualizeParticleSimulation2D(std::string fileName);
```

**Example:**
```cpp
// Data already saved by simulation code
Visualizer::VisualizeParticleSimulation2D("collision_sim_100_balls.txt");
```

**Data Format:** File should contain particle positions/velocities over time (format depends on specific visualizer).

#### 3D Particle Simulation

```cpp
static void VisualizeParticleSimulation3D(std::string fileName);
```

**Example:**
```cpp
// 3D N-body simulation
Visualizer::VisualizeParticleSimulation3D("nbody_galaxy.txt");
```

---

## Usage Examples

### Complete Workflow Examples

#### Example 1: Function Analysis

```cpp
#include "tools/Visualizer.h"
#include "base/RealFunction.h"

void analyzeFunctionVisually() {
    // Define function
    auto gaussian = [](Real x) {
        return std::exp(-x*x);
    };
    RealFunction f(gaussian);
    
    // Visualize over interval
    Visualizer::VisualizeRealFunction(
        f,
        "Gaussian exp(-x²)",
        -3.0, 3.0,
        200,
        "gaussian.txt"
    );
    
    // Compare with approximations
    auto taylor2 = [](Real x) { return 1 - x*x; };
    auto taylor4 = [](Real x) { return 1 - x*x + x*x*x*x/2.0; };
    
    RealFunction t2(taylor2), t4(taylor4);
    std::vector<IRealFunction*> funcs = {&f, &t2, &t4};
    std::vector<std::string> legend = {"exp(-x²)", "Taylor O(2)", "Taylor O(4)"};
    
    Visualizer::VisualizeMultiRealFunction(
        funcs,
        "Taylor Approximations",
        legend,
        -2.0, 2.0,
        150,
        "taylor_comparison.txt"
    );
}
```

#### Example 2: Vector Field Visualization

```cpp
#include "tools/Visualizer.h"
#include "base/VectorFunction.h"

void visualizeElectricField() {
    // Electric field of point charge at origin
    auto eField = [](const VectorN<Real, 2>& r) {
        Real x = r[0], y = r[1];
        Real rMag = std::sqrt(x*x + y*y);
        if (rMag < 1e-10) return VectorN<Real, 2>{0, 0};
        
        Real scale = 1.0 / (rMag * rMag * rMag);
        return VectorN<Real, 2>{x * scale, y * scale};
    };
    
    VectorFunction<2> field(eField);
    
    Visualizer::VisualizeVectorField2DCartesian(
        field,
        "Electric Field (Point Charge)",
        -5.0, 5.0, 25,
        -5.0, 5.0, 25,
        "electric_field.txt"
    );
}
```

#### Example 3: ODE Solution Phase Portrait

```cpp
#include "tools/Visualizer.h"
#include "algorithms/ODESystemSolver.h"

void visualizePendulum() {
    // Pendulum: θ'' + sin(θ) = 0
    // State: [θ, ω]
    auto pendulum = [](Real t, const VectorN<Real, 2>& state) {
        Real theta = state[0];
        Real omega = state[1];
        return VectorN<Real, 2>{omega, -std::sin(theta)};
    };
    
    ODESystem<2> system(pendulum);
    VectorN<Real, 2> initState{M_PI/4, 0.0};  // 45° initial angle
    
    ODESystemFixedStepSolver<2> solver(system);
    auto sol = solver.integrate_RK4(initState, 0.0, 10.0, 0.01);
    
    // Time series of angle
    Visualizer::VisualizeODESysSolCompAsFunc(
        sol, 0, "Angle vs Time", "theta_t.txt"
    );
    
    // Phase portrait (ω vs θ)
    Visualizer::VisualizeODESysSolAsParamCurve2(
        sol, 0, 1, "Phase Portrait", "phase.txt"
    );
    
    // Both variables
    std::vector<std::string> labels = {"θ", "ω"};
    Visualizer::VisualizeODESysSolAsMultiFunc(
        sol, "Pendulum State", labels, "state.txt"
    );
}
```

#### Example 4: 3D Curves

```cpp
void visualize3DCurves() {
    // Trefoil knot
    auto trefoil = [](Real t) {
        Real x = std::sin(t) + 2*std::sin(2*t);
        Real y = std::cos(t) - 2*std::cos(2*t);
        Real z = -std::sin(3*t);
        return VectorN<Real, 3>{x, y, z};
    };
    
    RealToVectorFunction<3> knot(trefoil);
    
    Visualizer::VisualizeParamCurve3D(
        knot,
        "Trefoil Knot",
        0.0, 2*M_PI,
        500,
        "trefoil.txt"
    );
}
```

---

## Data Export Formats

All visualization methods export data to text files in the `results/` directory. The format depends on the data type:

### Real Function Format

```
REAL_FUNCTION
<title>
x1: <start>
x2: <end>
NumPoints: <count>
<x1> <f(x1)>
<x2> <f(x2)>
...
```

### Multi-Function Format

```
MULTI_REAL_FUNCTION
<title>
NumFunctions: <count>
NumPoints: <points>
<legend1>
<legend2>
...
<x1> <f1(x1)> <f2(x1)> ...
<x2> <f1(x2)> <f2(x2)> ...
...
```

### Parametric Curve 2D Format

```
PARAMETRIC_CURVE_2D
<title>
t1: <start>
t2: <end>
NumPoints: <count>
<x1> <y1>
<x2> <y2>
...
```

### Parametric Curve 3D Format

```
PARAMETRIC_CURVE_3D
<title>
t1: <start>
t2: <end>
NumPoints: <count>
<x1> <y1> <z1>
<x2> <y2> <z2>
...
```

### Scalar Field Format

```
SCALAR_FUNCTION_2D
<title>
x1: <xmin> x2: <xmax> NumPointsX: <nx>
y1: <ymin> y2: <ymax> NumPointsY: <ny>
<x1> <y1> <f(x1,y1)>
<x1> <y2> <f(x1,y2)>
...
```

### Vector Field Format

```
VECTOR_FIELD_2D
<title>
x1: <xmin> x2: <xmax> NumPointsX: <nx>
y1: <ymin> y2: <ymax> NumPointsY: <ny>
<x1> <y1> <Fx(x1,y1)> <Fy(x1,y1)>
<x1> <y2> <Fx(x1,y2)> <Fy(x1,y2)>
...
```

**Note:** All formats use space-separated values with high precision (15 decimal places by default).

---

## External Visualizers

### Windows Visualizers (WPF)

Located in `tools/visualizers/win/WPF/`:

- **RealFuncVisualizer**: 2D function plots with zoom, pan, grid
- **SurfaceVisualizer**: 3D surface plots with rotation, lighting
- **ParametricCurve2DVisualizer**: 2D parametric curve plots
- **ParametricCurve3DVisualizer**: 3D curve plots with camera controls
- **VectorField2DVisualizer**: 2D arrow plots with color mapping
- **VectorField3DVisualizer**: 3D vector field visualization
- **Particle2DVisualizer**: 2D particle animation
- **Particle3DVisualizer**: 3D particle/N-body simulation

**Features:** Interactive zoom, pan, rotation, export to image, grid/axis controls.

### Linux Visualizers

Located in `tools/visualizers/linux/`:

- **FLTK-based**: Lightweight, fast rendering
- **Qt-based**: Full-featured with advanced controls

**Status:** Partial implementation. Some visualizers planned but not yet integrated.

### Custom Visualizers

You can create custom visualizers by:

1. Reading the exported text format
2. Parsing header information (title, ranges, dimensions)
3. Rendering data using your preferred graphics library

**Example frameworks:** OpenGL, SDL, Qt, FLTK, Python (matplotlib/plotly), JavaScript (D3.js/Plotly.js).

---

## Platform-Specific Behavior

### Windows

- All visualization methods fully supported
- WPF visualizers launched via `system()` call
- Paths resolved automatically

### Linux

- Most methods supported
- ODE phase portrait visualization (`VisualizeODESysSolAsParamCurve2/3`) prints "Not implemented" message
- Data still exported; manual visualization possible

### Cross-Platform Data Export

All `Visualizer` methods export data files regardless of platform. Even if the visualizer launch fails, you can:

1. Find data in `results/` directory
2. Use external tools (MATLAB, Python, gnuplot, etc.)
3. Implement custom visualizer for your platform

---

## Best Practices

### 1. Choose Appropriate Sampling Density

```cpp
// Too sparse - misses details
Visualizer::VisualizeRealFunction(f, "Undersampled", 0, 10, 10, "bad.txt");

// Good - captures features
Visualizer::VisualizeRealFunction(f, "Well Sampled", 0, 10, 200, "good.txt");

// Overkill - wastes resources
Visualizer::VisualizeRealFunction(f, "Oversampled", 0, 10, 10000, "slow.txt");
```

**Guidelines:**
- Smooth functions: 100-200 points
- Oscillatory functions: 10-20 points per period
- Discontinuities: Dense sampling near jumps

### 2. Use Descriptive Titles and Filenames

```cpp
// Poor
Visualizer::VisualizeRealFunction(f, "f", 0, 1, 100, "out.txt");

// Good
Visualizer::VisualizeRealFunction(
    f, 
    "Velocity Profile (Re=1000)", 
    0, 1, 100, 
    "velocity_profile_re1000.txt"
);
```

### 3. Separate Files for Different Visualizations

```cpp
// Organize results by category
Visualizer::VisualizeRealFunction(f, "...", 0, 1, 100, "analysis/case1_func.txt");
Visualizer::VisualizeScalarFunc2DCartesian(s, "...", -1, 1, 50, -1, 1, 50, "analysis/case1_field.txt");
```

### 4. Verify Data Before Visualization

```cpp
// Check ODE solution validity before visualizing
if (sol.IsValid()) {
    Visualizer::VisualizeODESysSolAsMultiFunc(sol, "...", legend, "ode_sol.txt");
} else {
    std::cerr << "ODE solution invalid - visualization skipped" << std::endl;
}
```

### 5. Batch Visualizations for Comparison

```cpp
// Compare numerical methods
std::vector<ODESystemSolution> solutions = {solEuler, solRK4, solRKCK};
std::vector<std::string> methods = {"Euler", "RK4", "RKCK"};

for (size_t i = 0; i < solutions.size(); ++i) {
    std::string filename = "method_" + methods[i] + ".txt";
    Visualizer::VisualizeODESysSolAsParamCurve2(
        solutions[i], 0, 1, methods[i] + " Phase Portrait", filename
    );
}
```

---

## Performance Considerations

### File I/O

- **Serialization time**: Proportional to number of points
- **File size**: ~30 bytes per point (with 15-digit precision)
- **Large datasets**: Consider reducing precision or sampling

### External Process Launch

```cpp
// system() call blocks until visualizer closes
Visualizer::VisualizeRealFunction(f, "...", 0, 1, 100, "data.txt");
// Execution continues only after user closes visualizer window
```

**Workaround for non-blocking:**
- Export data without launching visualizer
- Use `Serializer` directly for data-only export
- Launch visualizer separately/manually

### Memory Usage

```cpp
// High memory for dense sampling
int nx = 1000, ny = 1000;  // 1 million points
Visualizer::VisualizeScalarFunc2DCartesian(
    surface, "Dense Grid", -10, 10, nx, -10, 10, ny, "dense.txt"
);
// Requires ~30 MB file + evaluation overhead
```

**Recommendation:** Start with coarse sampling, refine as needed.

---

## Integration with Other Components

### With Serializer

```cpp
// Manual export for custom processing
Serializer::SaveRealFunc(f, "Function", 0, 10, 100, "results/data.txt");

// Then use Visualizer for same data (re-exports)
Visualizer::VisualizeRealFunction(f, "Function", 0, 10, 100, "data.txt");
```

**Tip:** Use `Serializer` directly when you only need data export without visualization.

### With ODE Solvers

```cpp
// Solve system
ODESystemFixedStepSolver<4> solver(system);
auto sol = solver.integrate_RK4(initState, 0.0, 100.0, 0.01);

// Visualize multiple aspects
Visualizer::VisualizeODESysSolAsMultiFunc(sol, "All Vars", labels, "all.txt");
Visualizer::VisualizeODESysSolAsParamCurve2(sol, 0, 1, "Phase", "phase.txt");
Visualizer::VisualizeODESysSolCompAsFunc(sol, 0, "X vs t", "x_t.txt");
```

### With Function Analysis

```cpp
#include "algorithms/FunctionsAnalyzer.h"

// Find and visualize extrema
auto criticalPoints = FunctionsAnalyzer::FindCriticalPoints(f, a, b, dx);

// Visualize function with critical points marked
Vector<Real> xCritical(criticalPoints.size());
for (size_t i = 0; i < criticalPoints.size(); ++i) {
    xCritical[i] = criticalPoints[i].position;
}

Visualizer::VisualizeRealFunction(f, "Function", a, b, 200, "func.txt");
Visualizer::VisualizeRealFunction(f, "Critical Points", xCritical, "critical.txt");
```

---

## Common Issues and Solutions

### Issue: Visualizer Doesn't Launch

**Possible causes:**
1. Path configuration incorrect
2. Visualizer executable missing/not compiled
3. Platform mismatch (Linux visualizers not available)

**Solutions:**
```cpp
// Check path configuration in MMLBase.h
// Verify GetRealFuncVisualizerPath() returns valid path

// Data is still exported - use external tool:
// Python: plt.loadtxt("results/data.txt"); plt.plot(...)
// gnuplot: plot "results/data.txt" using 1:2 with lines
// MATLAB: data = load("results/data.txt"); plot(data(:,1), data(:,2))
```

### Issue: File Not Found Error

**Cause:** Relative paths assume working directory contains `results/` folder.

**Solution:**
```cpp
// Ensure results directory exists
#include <filesystem>
std::filesystem::create_directory("results");

// Or use absolute paths in Serializer configuration
```

### Issue: Visualization Shows Nothing

**Cause:** Function domain outside visualizer viewport.

**Solution:**
- Use appropriate range for data
- Check visualizer zoom/pan controls
- Verify function is not NaN/Inf in range

### Issue: Poor Quality Curves

**Cause:** Insufficient sampling density.

**Solution:**
```cpp
// Increase numPoints
Visualizer::VisualizeRealFunction(f, "Better", 0, 10, 500, "better.txt");  // Was 50

// Or use adaptive sampling via Serializer with custom points
```

---

## See Also

- **[Serializer](Serializer.md)**: Direct data export without visualization
- **[ConsolePrinter](ConsolePrinter.md)**: Terminal-based output for debugging
- **[ODE Solvers](README_Algorithms.md#ode-solvers)**: Solve systems for visualization
- **[FunctionsAnalyzer](README_Algorithms.md#function-analysis)**: Analyze functions before visualizing

---

## Summary

The **Visualizer** class provides:
- ✅ **10+ visualization types**: Functions, fields, curves, surfaces, particles
- ✅ **Platform integration**: Windows WPF visualizers, Linux support planned
- ✅ **Data portability**: Standard text formats readable by any tool
- ✅ **Batch visualization**: Compare multiple objects simultaneously
- ✅ **ODE integration**: Direct visualization of ODE solutions
- ✅ **External tool compatibility**: Export data for Python/MATLAB/gnuplot

**Key Takeaway:** Visualizer bridges computational results and graphical display, making it easy to verify algorithms, present results, and debug simulations through visual inspection.
