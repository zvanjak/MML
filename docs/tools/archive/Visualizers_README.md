# Visualizers - Graphical Display of Mathematical Objects

**Header:** `mml/tools/Visualizer.h`  
**Namespace:** `MML`  
**Type:** Static utility class

---

## Overview

The Visualizer class provides a unified interface for creating graphical visualizations of mathematical objects. It automatically serializes data using the Serializer class and launches external plotting applications to display the results. All methods are static and platform-aware.

### Key Features

✅ **Real Functions** - Plot 1D functions and comparisons  
✅ **Multi-Functions** - Overlay multiple functions on one plot  
✅ **Scalar Fields** - 3D surface plots of 2D scalar functions  
✅ **Vector Fields** - 2D and 3D arrow field visualizations  
✅ **Parametric Curves** - 2D and 3D trajectory plotting  
✅ **ODE Solutions** - Time-series and phase portraits  
✅ **Particle Simulations** - 2D and 3D animations  
✅ **Automatic Serialization** - No manual file management  
✅ **External Viewers** - Launches platform-appropriate plotting tools  

---

## Quick Start

### Visualize a Simple Function

```cpp
#include "tools/Visualizer.h"
using namespace MML;

// Define a function
auto f = [](Real x) { return std::sin(x); };
RealFunction sine(f);

// Visualize immediately
Visualizer::VisualizeRealFunction(sine, "sin(x)", 
                                  0.0, 2*M_PI, 200, 
                                  "sine.txt");
// Window opens automatically showing the plot
```

### Compare Multiple Functions

```cpp
auto sine = [](Real x) { return std::sin(x); };
auto cosine = [](Real x) { return std::cos(x); };
RealFunction f1(sine), f2(cosine);

std::vector<IRealFunction*> funcs = {&f1, &f2};
std::vector<std::string> legend = {"sin(x)", "cos(x)"};

Visualizer::VisualizeMultiRealFunction(
    funcs, "Trigonometric Functions", legend,
    0, 2*M_PI, 200, "trig.txt");
```

### Visualize 3D Parametric Curve

```cpp
// Helix curve
auto helix = [](Real t) -> Vector<Real> {
    return Vector<Real>{std::cos(t), std::sin(t), 0.2*t};
};
VectorFunction<3> curve(helix);

Visualizer::VisualizeParamCurve3D(
    curve, "Helix", 0, 4*M_PI, 300, "helix.txt");
```

---

## Architecture

### How Visualizer Works

1. **Serialize**: Calls `Serializer::Save*()` to write data to file
2. **Locate Viewer**: Uses platform paths to find external plotting application
3. **Launch**: Executes viewer with file path as argument
4. **Display**: External viewer opens window and renders plot

### Path Configuration

Visualizer uses getter functions for cross-platform path resolution:

```cpp
_pathResultFiles          // Where data files are saved
_pathRealFuncViz          // 1D function plotter
_pathSurfaceViz           // 3D surface plotter
_pathParametricCurve2DViz // 2D curve plotter
_pathParametricCurve3DViz // 3D curve plotter
_pathVectorField2DViz     // 2D vector field plotter
_pathVectorField3DViz     // 3D vector field plotter
_pathParticle2DViz        // 2D particle animation viewer
_pathParticle3DViz        // 3D particle animation viewer
```

These paths are configured in `MMLBase.h` and adapt to the platform.

---

## API Reference

### Real Function Visualization

#### VisualizeRealFunction (Equally-Spaced)

Plot a real function over a uniform grid.

```cpp
static void VisualizeRealFunction(
    const IRealFunction& f,    // Function to plot
    std::string title,         // Plot title
    Real x1, Real x2,          // Domain [x1, x2]
    int numPoints,             // Number of sample points
    std::string fileName       // Output filename (data + plot)
);
```

**Example:**
```cpp
auto cubic = [](Real x) { return x*x*x - 3*x + 1; };
RealFunction f(cubic);

Visualizer::VisualizeRealFunction(
    f, "Cubic Function", -2, 2, 100, "cubic.txt");
```

#### VisualizeRealFunction (Arbitrary Points)

Plot a function at specific x-coordinates.

```cpp
static void VisualizeRealFunction(
    const IRealFunction& f,
    std::string title,
    Vector<Real> points,       // X-coordinates to evaluate
    std::string fileName
);
```

**Example (Adaptive Sampling):**
```cpp
// More points near x=0 where function changes rapidly
Vector<Real> points(100);
for (int i = 0; i < 100; i++) {
    Real t = i / 99.0;
    points[i] = std::pow(t, 2) * 4 - 2;  // Quadratic spacing
}

Visualizer::VisualizeRealFunction(
    f, "Adaptive Sampling", points, "adaptive.txt");
```

---

### Multi-Function Visualization

#### VisualizeMultiRealFunction (Function Pointers)

Plot multiple functions on the same axes.

```cpp
static void VisualizeMultiRealFunction(
    std::vector<IRealFunction*> funcs,    // Functions to plot
    std::string title,
    std::vector<std::string> func_legend, // Legend labels
    Real x1, Real x2,
    int numPoints,
    std::string fileName
);
```

**Example (Polynomial Basis):**
```cpp
auto p0 = [](Real x) { return 1.0; };
auto p1 = [](Real x) { return x; };
auto p2 = [](Real x) { return x*x; };
auto p3 = [](Real x) { return x*x*x; };

RealFunction f0(p0), f1(p1), f2(p2), f3(p3);
std::vector<IRealFunction*> funcs = {&f0, &f1, &f2, &f3};
std::vector<std::string> legend = {"1", "x", "x²", "x³"};

Visualizer::VisualizeMultiRealFunction(
    funcs, "Polynomial Basis", legend, 
    -1, 1, 100, "polynomials.txt");
```

#### VisualizeMultiRealFunction (Interpolated Functions)

Specialized overloads for interpolated function types:

```cpp
// Linear interpolation
static void VisualizeMultiRealFunction(
    std::vector<LinearInterpRealFunc> funcs, ...);

// Polynomial interpolation
static void VisualizeMultiRealFunction(
    std::vector<PolynomInterpRealFunc> funcs, ...);

// Spline interpolation
static void VisualizeMultiRealFunction(
    std::vector<SplineInterpRealFunc> funcs, ...);
```

#### VisualizeMultiRealFunctionSeparately

Visualize multiple functions with different domains.

```cpp
static void VisualizeMultiRealFunctionSeparately(
    std::vector<LinearInterpRealFunc> funcs,
    std::string title,
    std::vector<std::string> func_legend,
    Real x1, Real x2,
    int numPoints,
    std::string fileName
);
```

**Use Case:** When functions have different valid ranges but you want to compare them.

---

### Scalar Field Visualization

#### VisualizeScalarFunc2DCartesian

Display a 2D scalar field as a 3D surface.

```cpp
static void VisualizeScalarFunc2DCartesian(
    const IScalarFunction<2>& func,
    std::string title,
    Real x1, Real x2, int numPointsX,  // X-axis grid
    Real y1, Real y2, int numPointsY,  // Y-axis grid
    std::string fileName
);
```

**Example (Gaussian Peak):**
```cpp
auto gaussian = [](const Vector<Real>& v) {
    Real x = v[0], y = v[1];
    return std::exp(-(x*x + y*y));
};
ScalarFunction<2> f(gaussian);

Visualizer::VisualizeScalarFunc2DCartesian(
    f, "2D Gaussian", 
    -3, 3, 50,   // 50 points in x
    -3, 3, 50,   // 50 points in y
    "gaussian2d.txt");
```

---

### Vector Field Visualization

#### VisualizeVectorField2DCartesian

Display 2D vector field with arrows.

```cpp
static void VisualizeVectorField2DCartesian(
    const IVectorFunction<2>& func,
    std::string title,
    Real x1, Real x2, int numPointsX,
    Real y1, Real y2, int numPointsY,
    std::string fileName
);
```

**Example (Vortex):**
```cpp
auto vortex = [](const Vector<Real>& v) -> Vector<Real> {
    Real x = v[0], y = v[1];
    return Vector<Real>{-y, x};  // Counterclockwise rotation
};
VectorFunction<2> field(vortex);

Visualizer::VisualizeVectorField2DCartesian(
    field, "Vortex Field", 
    -2, 2, 20, -2, 2, 20, "vortex.txt");
```

#### VisualizeVectorField3DCartesian

Display 3D vector field.

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

---

### Parametric Curve Visualization

#### VisualizeParamCurve2D

Display a 2D parametric curve.

```cpp
static void VisualizeParamCurve2D(
    const IRealToVectorFunction<2>& f,
    std::string title,
    Real t1, Real t2,          // Parameter range
    int numPoints,             // Resolution
    std::string fileName
);
```

**Example (Lissajous Figure):**
```cpp
auto lissajous = [](Real t) -> Vector<Real> {
    return Vector<Real>{std::sin(3*t), std::sin(2*t)};
};
VectorFunction<2> curve(lissajous);

Visualizer::VisualizeParamCurve2D(
    curve, "Lissajous 3:2", 0, 2*M_PI, 500, "lissajous.txt");
```

#### VisualizeMultiParamCurve2D (Direct)

Display multiple 2D curves on same plot.

```cpp
static void VisualizeMultiParamCurve2D(
    std::vector<IRealToVectorFunction<2>*> curves,
    std::string title,
    Real t1, Real t2,
    int numPoints,
    std::string fileName
);
```

**Example (Family of Circles):**
```cpp
std::vector<IRealToVectorFunction<2>*> circles;
std::vector<std::string> legend;

for (int r = 1; r <= 3; r++) {
    auto circle = new VectorFunction<2>([r](Real t) -> Vector<Real> {
        return Vector<Real>{r*std::cos(t), r*std::sin(t)};
    });
    circles.push_back(circle);
    legend.push_back("r=" + std::to_string(r));
}

Visualizer::VisualizeMultiParamCurve2D(
    circles, "Concentric Circles", legend,
    0, 2*M_PI, 100, "circles.txt");
```

#### VisualizeMultiParamCurve2D (Pre-Serialized)

Display multiple curves from existing data files.

```cpp
static void VisualizeMultiParamCurve2D(
    std::vector<std::string> fileNames  // Already serialized data
);
```

**Use Case:** When you have data from previous runs or external sources.

#### VisualizeParamCurve3D

Display a 3D parametric curve.

```cpp
static void VisualizeParamCurve3D(
    const IRealToVectorFunction<3>& f,
    std::string title,
    Real t1, Real t2,
    int numPoints,
    std::string fileName
);
```

**Example (Trefoil Knot):**
```cpp
auto trefoil = [](Real t) -> Vector<Real> {
    return Vector<Real>{
        std::sin(t) + 2*std::sin(2*t),
        std::cos(t) - 2*std::cos(2*t),
        -std::sin(3*t)
    };
};
VectorFunction<3> knot(trefoil);

Visualizer::VisualizeParamCurve3D(
    knot, "Trefoil Knot", 0, 2*M_PI, 500, "trefoil.txt");
```

#### VisualizeMultiParamCurve3D

Display multiple 3D curves.

```cpp
static void VisualizeMultiParamCurve3D(
    std::vector<IRealToVectorFunction<3>*> curves,
    std::string title,
    Real t1, Real t2,
    int numPoints,
    std::string fileName
);
```

---

### ODE Solution Visualization

#### VisualizeODESysSolCompAsFunc

Plot a single component vs time.

```cpp
static void VisualizeODESysSolCompAsFunc(
    const ODESystemSolution& sol,
    int compInd,               // Component index
    std::string title,
    std::string fileName
);
```

**Example (Pendulum Angle):**
```cpp
// After solving pendulum ODE
ODESystemSolution solution = solver.solve(pendulum, init_cond);

// Plot angle vs time (component 0)
Visualizer::VisualizeODESysSolCompAsFunc(
    solution, 0, "Pendulum Angle θ(t)", "angle.txt");
```

#### VisualizeODESysSolAsMultiFunc

Plot all components vs time.

```cpp
static void VisualizeODESysSolAsMultiFunc(
    const ODESystemSolution& sol,
    std::string title,
    std::vector<std::string> legend,   // Component names
    std::string fileName
);
```

**Example (State Variables):**
```cpp
std::vector<std::string> legend = {"x", "y", "vx", "vy"};
Visualizer::VisualizeODESysSolAsMultiFunc(
    solution, "2D Projectile Motion", legend, "projectile.txt");
```

#### VisualizeODESysSolAsParamCurve2

Plot two components against each other (phase portrait).

```cpp
static void VisualizeODESysSolAsParamCurve2(
    const ODESystemSolution& sol,
    int ind1, int ind2,        // Component indices
    std::string title,
    std::string fileName
);
```

**Example (Phase Space):**
```cpp
// Plot position vs velocity
Visualizer::VisualizeODESysSolAsParamCurve2(
    solution, 0, 1, "Phase Portrait: x vs v", "phase.txt");
```

#### VisualizeODESysSolAsParamCurve3

Plot three components as 3D trajectory.

```cpp
static void VisualizeODESysSolAsParamCurve3(
    const ODESystemSolution& sol,
    int ind1, int ind2, int ind3,
    std::string title,
    std::string fileName
);
```

**Example (Lorenz Attractor):**
```cpp
Visualizer::VisualizeODESysSolAsParamCurve3(
    lorenz_solution, 0, 1, 2,
    "Lorenz Attractor", "lorenz.txt");
```

---

### Particle Simulation Visualization

#### VisualizeParticleSimulation2D

Display 2D particle animation.

```cpp
static void VisualizeParticleSimulation2D(
    std::string fileName  // Pre-serialized animation data
);
```

**Workflow:**
```cpp
// 1. Simulate particles
std::vector<std::vector<Vector2Cartesian>> trajectory;
// ... run simulation ...

// 2. Serialize to file
Serializer::SaveParticleSimulation2D(
    "animation.txt", numBalls, width, height, trajectory);

// 3. Visualize
Visualizer::VisualizeParticleSimulation2D("animation.txt");
```

#### VisualizeParticleSimulation3D

Display 3D particle animation.

```cpp
static void VisualizeParticleSimulation3D(
    std::string fileName
);
```

---

## Usage Examples

### Example 1: Function Comparison Study

```cpp
// Compare numerical integration methods
auto exact = [](Real x) { 
    return std::exp(-x*x); 
};

// Generate approximate solution
Vector<Real> x_data = {0, 0.5, 1.0, 1.5, 2.0};
Vector<Real> y_data(5);
for (int i = 0; i < 5; i++)
    y_data[i] = exact(x_data[i]);

RealFunction exact_func(exact);
LinearInterpRealFunc linear(x_data, y_data);
SplineInterpRealFunc spline(x_data, y_data);

std::vector<IRealFunction*> funcs = {&exact_func};
std::vector<std::string> legend = {"Exact"};

// Plot exact solution
Visualizer::VisualizeMultiRealFunction(
    funcs, "Gaussian Function", legend, 
    0, 2, 200, "gaussian_exact.txt");

// Compare interpolation methods
std::vector<LinearInterpRealFunc> linear_vec = {linear};
Visualizer::VisualizeMultiRealFunction(
    linear_vec, "Linear Interpolation", {"Linear"},
    0, 2, 200, "linear_interp.txt");

std::vector<SplineInterpRealFunc> spline_vec = {spline};
Visualizer::VisualizeMultiRealFunction(
    spline_vec, "Spline Interpolation", {"Spline"},
    0, 2, 200, "spline_interp.txt");
```

### Example 2: ODE Solution Analysis

```cpp
// Define ODE system (damped oscillator)
auto damped_oscillator = [](Real t, const Vector<Real>& y) -> Vector<Real> {
    Real x = y[0], v = y[1];
    Real omega = 2.0, gamma = 0.1;
    return Vector<Real>{v, -omega*omega*x - 2*gamma*v};
};

ODESystem system(damped_oscillator, 2);
Vector<Real> init{1.0, 0.0};  // x(0)=1, v(0)=0

// Solve
ODESolver solver(/*...*/);
ODESystemSolution sol = solver.solve(system, init, 0, 20, 0.01);

// Visualize position vs time
Visualizer::VisualizeODESysSolCompAsFunc(
    sol, 0, "Position x(t)", "position.txt");

// Visualize velocity vs time
Visualizer::VisualizeODESysSolCompAsFunc(
    sol, 1, "Velocity v(t)", "velocity.txt");

// Visualize both on same plot
std::vector<std::string> legend = {"Position", "Velocity"};
Visualizer::VisualizeODESysSolAsMultiFunc(
    sol, "Damped Oscillator", legend, "both.txt");

// Phase portrait
Visualizer::VisualizeODESysSolAsParamCurve2(
    sol, 0, 1, "Phase Portrait: x vs v", "phase.txt");
```

### Example 3: 3D Trajectory Visualization

```cpp
// Charged particle in magnetic field
auto magnetic_field = [](Real t, const Vector<Real>& y) -> Vector<Real> {
    Real x = y[0], y_pos = y[1], z = y[2];
    Real vx = y[3], vy = y[4], vz = y[5];
    
    // Uniform magnetic field in z-direction
    Real B = 1.0, q_over_m = 1.0;
    Real ax = q_over_m * vy * B;
    Real ay = -q_over_m * vx * B;
    Real az = 0;
    
    return Vector<Real>{vx, vy, vz, ax, ay, az};
};

ODESystem system(magnetic_field, 6);
Vector<Real> init{0, 0, 0, 1, 0, 1};  // Initial position and velocity

ODESystemSolution sol = solver.solve(system, init, 0, 10, 0.01);

// Visualize 3D helical trajectory
Visualizer::VisualizeODESysSolAsParamCurve3(
    sol, 0, 1, 2,  // x, y, z components
    "Charged Particle Trajectory", "helix.txt");
```

### Example 4: Vector Field and Trajectories

```cpp
// Define vector field
auto flow_field = [](const Vector<Real>& v) -> Vector<Real> {
    Real x = v[0], y = v[1];
    return Vector<Real>{-y, x};  // Rotation field
};
VectorFunction<2> field(flow_field);

// Visualize the field
Visualizer::VisualizeVectorField2DCartesian(
    field, "Rotation Field", 
    -2, 2, 15, -2, 2, 15, "field.txt");

// Now simulate particles in this field and visualize trajectories
// ... (create trajectories) ...

// Visualize multiple particle paths
std::vector<IRealToVectorFunction<2>*> paths;
// ... (add paths) ...

Visualizer::VisualizeMultiParamCurve2D(
    paths, "Particle Trajectories in Flow Field",
    0, 10, 200, "trajectories.txt");
```

---

## Platform Notes

### Windows
All visualization methods fully supported. Uses external viewers configured in `MMLBase.h`.

### Linux/macOS
Most methods supported. Some methods check platform with `#ifdef _WIN32` and print "Not implemented" message on other platforms.

**Platform-Specific Methods:**
```cpp
VisualizeODESysSolAsParamCurve2  // Windows only
VisualizeODESysSolAsParamCurve3  // Windows only
```

---

## Best Practices

### 1. Choose Appropriate Resolution

```cpp
// For smooth curves
Visualizer::VisualizeRealFunction(f, title, x1, x2, 500, file);

// For quick preview
Visualizer::VisualizeRealFunction(f, title, x1, x2, 50, file);

// For very high precision (slow!)
Visualizer::VisualizeRealFunction(f, title, x1, x2, 5000, file);
```

### 2. Use Descriptive Titles

```cpp
// Good
Visualizer::VisualizeRealFunction(
    f, "sin(x)", 0, 2*M_PI, 200, "sine.txt");

// Better
Visualizer::VisualizeRealFunction(
    f, "Sine function: period=2π, amplitude=1", 
    0, 2*M_PI, 200, "sine.txt");
```

### 3. Organize Output Files

```cpp
// Group by category
Visualizer::VisualizeRealFunction(f, title, x1, x2, n, "functions/sine.txt");
Visualizer::VisualizeODESysSolAsMultiFunc(sol, title, legend, "ode/pendulum.txt");
Visualizer::VisualizeVectorField2DCartesian(field, title, x1,x2,nx, y1,y2,ny, "fields/gradient.txt");
```

### 4. Balance Grid Resolution

```cpp
// 2D scalar field: 50×50 = 2500 points (good balance)
Visualizer::VisualizeScalarFunc2DCartesian(
    f, title, -3, 3, 50, -3, 3, 50, "field.txt");

// 3D vector field: 20×20×20 = 8000 points (reasonable)
Visualizer::VisualizeVectorField3DCartesian(
    field, title, -2, 2, 20, -2, 2, 20, -2, 2, 20, "field3d.txt");
```

### 5. Check Platform Before Specialized Features

```cpp
#ifdef _WIN32
    Visualizer::VisualizeODESysSolAsParamCurve3(
        sol, 0, 1, 2, title, "trajectory.txt");
#else
    // Fallback: save data only
    Serializer::SaveODESolAsParametricCurve3D(
        sol, "trajectory.txt", 0, 1, 2, title);
    std::cout << "Data saved. Use external viewer." << std::endl;
#endif
```

---

## Troubleshooting

### Visualization Window Doesn't Open

**Possible Causes:**
1. External viewer not installed/configured
2. Incorrect paths in `MMLBase.h`
3. File write failed (check Serializer error)

**Solution:**
```cpp
// Check if serialization succeeded
auto result = Serializer::SaveRealFunc(f, title, x1, x2, n, file);
if (!result.success) {
    std::cerr << "Serialization failed: " << result.message << std::endl;
}

// Manually check if viewer exists
// Update paths in MMLBase.h if needed
```

### Plot Looks Jagged/Rough

**Cause:** Too few sample points

**Solution:**
```cpp
// Increase numPoints
Visualizer::VisualizeRealFunction(
    f, title, x1, x2, 1000, file);  // Was 100
```

### Very Slow Visualization

**Cause:** Too many points in 2D/3D grids

**Solution:**
```cpp
// Reduce grid resolution
// From 100×100 = 10,000 points
Visualizer::VisualizeScalarFunc2DCartesian(
    f, title, -5, 5, 100, -5, 5, 100, file);

// To 50×50 = 2,500 points (4x faster)
Visualizer::VisualizeScalarFunc2DCartesian(
    f, title, -5, 5, 50, -5, 5, 50, file);
```

### Memory Issues with Particle Simulations

**Cause:** Too many frames or particles

**Solution:**
```cpp
// Reduce frame count or subsample
int totalFrames = 10000;
int subsample = 10;  // Keep every 10th frame

std::vector<std::vector<Vector2Cartesian>> trajectory;
for (int i = 0; i < totalFrames; i += subsample) {
    trajectory.push_back(positions[i]);
}

Serializer::SaveParticleSimulation2D(file, numBalls, w, h, trajectory);
```

---

## Integration with Other MML Tools

### Serializer + Visualizer Workflow

```cpp
// Manual serialization (for debugging or custom formats)
auto result = Serializer::SaveRealFunc(f, title, x1, x2, n, file);
if (result.success) {
    Visualizer::VisualizeRealFunction(f, title, x1, x2, n, file);
}

// Or use Visualizer directly (it calls Serializer internally)
Visualizer::VisualizeRealFunction(f, title, x1, x2, n, file);
```

### ConsolePrinter + Visualizer

```cpp
// Print data table to console
TablePrinter<Real> printer;
printer.AddColumn("x", points);
printer.AddColumn("f(x)", values);
printer.Print();

// Also visualize graphically
Visualizer::VisualizeRealFunction(f, title, x1, x2, n, file);
```

---

## Advanced Topics

### Custom Viewer Integration

To add a new visualization backend:

1. Add path getter in `MMLBase.h`:
```cpp
inline std::string GetCustomVisualizerPath() {
    return "/path/to/custom/viewer";
}
```

2. Add method in `Visualizer` class:
```cpp
static void VisualizeCustom(/*...*/) {
    std::string name = _pathResultFiles + fileName;
    Serializer::SaveCustomFormat(/*...*/);
    
    std::string command = GetCustomVisualizerPath() + " \"" + name + "\"";
    system(command.c_str());
}
```

### Batch Visualization

```cpp
std::vector<std::string> files;

// Generate multiple visualizations
for (int i = 0; i < 10; i++) {
    Real freq = i + 1;
    auto f = [freq](Real x) { return std::sin(freq * x); };
    RealFunction func(f);
    
    std::string filename = "sine_" + std::to_string(i) + ".txt";
    Visualizer::VisualizeRealFunction(
        func, "sin(" + std::to_string(freq) + "x)", 
        0, 2*M_PI, 200, filename);
    
    files.push_back(filename);
}
```

---

## See Also

- **Serializer** - Underlying data persistence layer
- **ConsolePrinter** - Console-based data display
- **ODESolver** - Generate solutions to visualize
- **Function Analysis** - Analyze functions before plotting

---

## Version History

**Current Version**
- Real function, multi-function, scalar field visualization
- 2D and 3D vector fields
- Parametric curves (2D and 3D)
- ODE solution visualization (time-series, phase portraits)
- Particle simulation animations
- Platform-aware viewer selection
- Automatic serialization with error handling
