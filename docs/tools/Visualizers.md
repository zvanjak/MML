# Visualizer - Data Visualization and Export

## Overview

The `Visualizer` class provides a high-level interface for visualizing mathematical objects and simulation results by exporting data to files and launching external visualization applications. It acts as a bridge between MML's computational components and platform-specific visualization tools.

## Header Organization

The public entry point remains `mml/tools/Visualizer.h`; existing `Visualizer::Visualize*`
calls and configuration functions are unchanged. It is now a small umbrella over four headers
under `mml/tools/visualizers/`:

| Header | Responsibility |
|--------|----------------|
| `VisualizerConfig.h` | Backend/root selection, executable catalog and executable paths |
| `VisualizerFiles.h` | Input lookup, sanitized output paths, containment and filesystem diagnostics |
| `VisualizerProcess.h` | Shell-free process execution, exit status and timeout handling |
| `VisualizerFacade.h` | Public mathematical/simulation export-and-launch workflows |

Configuration-only clients can include `mml/tools/visualizers/VisualizerConfig.h` without
loading serializers or Windows/POSIX process headers. The former `mml/MMLVisualizators.h`
path remains a compatibility include of that same configuration implementation. Both include
paths use the same runtime state; there is no second launcher or duplicated configuration.

The internal executable catalog is shared by availability checks, path getters and launch
workflows. Serialization remains in the existing serializer modules. The facade shares
save-and-launch and existing-file lookup logic rather than duplicating those workflows.

---

## Gallery

### Windows (WPF)

| Real Functions | Scalar Functions (Surfaces) |
|:--------------:|:---------------------------:|
| ![Real Function](../images/_Win%20-%20Real%20func%201.png) | ![Surface](../images/_Win%20-%20Scalar%20func%201.png) |
| ![Multi Function](../images/_Win%20-%20Real%20func%20-%20multi%201.png) | ![Surface 2](../images/_Win%20-%20Scalar%20func%202.png) |

| Parametric Curves | Vector Fields |
|:-----------------:|:-------------:|
| ![2D Curve](../images/_Win%20-%20Parametric%20curve%202D%201.png) | ![Vector Field 2D](../images/_Win%20-%20Vector%20field%202D%201.png) |
| ![3D Curve](../images/_Win%20-%20Parametric%20curve%203D%201.png) | ![Vector Field 2D 2](../images/_Win%20-%20Vector%20field%202D%202.png) |

| Particle Simulations |
|:--------------------:|
| ![Particle 2D](../images/_Win%20-%20Particle%20visualizer%202D.png) |
| ![Particle 3D](../images/_Win%20-%20Particle%20visualizer%203D.png) |

### Linux (Qt 6.9 + FLTK)

| Real Functions | Surfaces | Parametric Curves |
|:--------------:|:--------:|:-----------------:|
| ![Linux Real](../images/_Linux%20-%20Real%20function%20visualization.png) | ![Linux Surface](../images/_Linux%20-%20Surface%20visualization.png) | ![Linux 2D Curve](../images/_Linux%20-%20Parametric%20curve%202D.png) |

| 3D Curves | Vector Fields | Particle Simulations |
|:---------:|:-------------:|:--------------------:|
| ![Linux 3D Curve](../images/_Linux%20-%20Parametric%20curve%203D%20-%201.png) | ![Linux Vector](../images/_Linux%20-%20Vector%20field%202D.png) | ![Linux Particle](../images/_Linux%20-%20Particle%20visualizer%202D.png) |

---

## Platform Requirements

> **Platform-Specific Visualizers:**
> 
> | Platform | Framework | Use Case | Requirements |
> |----------|-----------|----------|-------------|
> | **Windows** | WPF | All visualizations | .NET 8 Runtime |
> | **Linux** | Qt 6.9 | All visualizations | Qt 6.9+ runtime |
> | **Linux** | FLTK | 2D only (lightweight) | FLTK 1.3+ |
> 
> **Windows:** WPF visualizers require .NET 8 Runtime (usually pre-installed on Windows 10/11, or download from [dotnet.microsoft.com](https://dotnet.microsoft.com/download)). Clone-mode binaries are bundled in `tools/visualizers/`; their source is maintained separately in the visualizer repository.
> 
> **Linux:** Two options available:
> - **Qt 6.9**: Full-featured visualizers for all types (2D/3D functions, curves, fields, particles)
> - **FLTK**: Lightweight alternative for 2D visualizations only
> 
> Source code: [MML_Linux_Visualizers](https://github.com/zvanjak/MML_Linux_Visualizers)

---

## Backend Selection API

MML provides a unified cross-platform backend selection system that allows you to choose which visualization framework to use at runtime.

### VisualizerBackend Enum

```cpp
enum class VisualizerBackend {
  WPF,   // Windows Presentation Foundation (Windows only, highest quality)
  Qt,    // Qt5/Qt6 (Cross-platform, feature-rich)
  FLTK,  // Fast Light Toolkit (Cross-platform, lightweight)
  Auto   // Platform-specific default (WPF on Windows, Qt on Linux)
};
```

### Backend Priority Order

The backend is determined in this order:
1. **Environment variable** `MML_VISUALIZER_BACKEND` (highest priority)
2. **Runtime API** setting via `SetVisualizerBackend()`
3. **Platform default** (WPF on Windows, Qt on Linux)

### Environment Variable

Override the backend without recompiling:

```bash
# Windows (Command Prompt)
set MML_VISUALIZER_BACKEND=Qt

# Windows (PowerShell)
$env:MML_VISUALIZER_BACKEND = "Qt"

# Linux
export MML_VISUALIZER_BACKEND=FLTK
```

Valid values: `WPF`, `Qt`, `QT5`, `QT6`, `FLTK`, `AUTO` (case-insensitive)

### Runtime API

```cpp
#include <mml/tools/visualizers/VisualizerConfig.h>  // or "MML.h" for single-header

using namespace MML;

// Check if a backend is available on this platform
bool available = IsBackendAvailable(VisualizerBackend::Qt);

// Set the backend at runtime
SetVisualizerBackend(VisualizerBackend::Qt);

// Get the currently active backend (never returns Auto, always resolves)
VisualizerBackend current = GetVisualizerBackend();
```

### Platform Defaults

| Platform | Default Backend | Available Backends |
|----------|-----------------|-------------------|
| Windows  | WPF | WPF, Qt, FLTK |
| Linux    | Qt  | Qt, FLTK |
| Other    | FLTK | FLTK |

### Example: Force Qt on Windows

```cpp
#include <MML.h>
using namespace MML;

int main() {
    // Use Qt instead of WPF on Windows for consistency with Linux
    if (IsBackendAvailable(VisualizerBackend::Qt)) {
        SetVisualizerBackend(VisualizerBackend::Qt);
    }
    
    // All subsequent visualizations use Qt
    RealFunction f([](Real x) { return std::sin(x); });
    Visualizer::VisualizeRealFunction(f, "sin(x)", 0, 2*M_PI, 100, "sine.mml");
    
    return 0;
}
```

### Notes

- If a requested backend is not available, `SetVisualizerBackend()` prints a warning and falls back to `Auto`.
- `GetVisualizerBackend()` always returns a concrete backend (never `Auto`).
- Backend availability requires at least one supported regular executable file
  in the configured bundle. Empty backend directories do not count. Partial
  installations are supported; each launch checks its exact executable.
- See `src/tryouts/test_linux_viz/test_backend_switching.cpp` for more examples.

---

**Key Features:**
- **Multiple visualization types**: Real functions, scalar/vector fields, parametric curves, surfaces, particles, ODE solutions
- **External visualizer integration**: Launches WPF (Windows) or Qt/FLTK (Linux) applications
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

## Highlights
- High-level bridge from MML computations to platform visualizers.
- Supports functions, fields, curves, surfaces, particles, and ODE solutions.
- Uses Serializer for robust data export; launches OS-specific viewers.
- Cross-platform path management for Windows/Linux visualizers.
- Multi-object visualization workflows for comparison and analysis.

---

## Working Examples

For complete, working code examples demonstrating all Visualizer capabilities, see:

📁 **[src/docs_demos/tools/docs_demo_visualizers.cpp](../../src/docs_demos/tools/docs_demo_visualizers.cpp)**

This file contains verified, compilable examples including:
- **Function Analysis**: Gaussian functions with Taylor approximations
- **Vector Fields**: Electric field visualization (2D and 3D)
- **ODE Solutions**: Pendulum time series, phase portraits, and multi-function plots
- **Parametric Curves**: Trefoil knot (3D), Lissajous figures (2D)
- **Scalar Functions**: 3D surface visualization
- **Lorenz Attractor**: Classic chaotic system as 3D trajectory

All examples are verified to compile and run correctly.

---

## Core Concepts

### Visualization Workflow

1. **Serialize**: Convert mathematical object to text file using `Serializer`
2. **Launch**: Execute external visualizer application with file path
3. **Display**: Visualizer app reads data and renders graphics

```cpp
// Example: Visualize a real function
IRealFunction& f = /* ... */;
Visualizer::VisualizeRealFunction(f, "sin(x)", 0.0, 2*M_PI, 100, "sine.mml");
// Creates: results/sine.txt
// Launches: visualizers/RealFuncVisualizer.exe "results/sine.mml"
```

### Path Configuration

Configure the **bundle root**, not the MML project root. The bundle contains
`win/`, `linux/`, and/or `mac/` with the existing backend layout beneath them.
Root selection is:

1. `SetVisualizerRoot()` runtime override.
2. `MML_VISUALIZER_ROOT` environment variable.
3. `visualizers/` under the current working directory.

No parent-directory search or source-checkout discovery is performed.
`MML_PROJECT_PATH` remains a result-file placement option; it does not select the
visualizer bundle. Existing clone-mode users should explicitly point the new root
to their `tools/visualizers` directory.

```cpp
const auto bundle = std::filesystem::current_path() / "visualizers";
SetVisualizerRoot(bundle);
const bool available = IsBackendAvailable(VisualizerBackend::Auto);
std::cout << "Visualizer bundle: " << GetVisualizerRoot() << '\n'
          << "Selected backend installed: " << std::boolalpha << available << '\n'
          << "Real-function executable: " << GetRealFuncVisualizerPath() << '\n';
```

Runnable counterpart: `Docs_Demo_VisualizerConfiguration()` in
`src/docs_demos/tools/docs_demo_visualizers.cpp`, selected with
`MML_DocsApp --visualizer-configuration`. This demo checks configuration without
opening a GUI.

Relative runtime roots are resolved when set; an empty root clears that override.
Environment and default roots are resolved when queried. Configure root/backend
before launching from other threads; concurrent mutation is not supported.
Path getters live in `mml/tools/visualizers/VisualizerConfig.h` and select Windows,
Linux, or macOS layout automatically. The legacy configuration include still exposes them.

### Launch results and literal paths

The development `MML_VisualizationApp` examples consume `IOResult` directly.
Failed exports/launches are reported to standard error with the operation name,
error category and message, plus the filename, nonzero exit code and timeout
state when present. Examples may continue to the next demonstration after a
reported failure; completion of the example list is not a guarantee that every
visualizer succeeded. Successful launches do not produce error diagnostics.
The `--help` path can be used to verify the CLI without launching a GUI.

Launchers execute the selected program directly, without a shell. Existing-file
paths containing spaces, parentheses, apostrophes, ampersands, dollar signs, or
other shell metacharacters are passed literally. Windows arguments are escaped
according to its command-line argument rules; POSIX arguments are passed through
`argv`. Embedded null characters cannot be represented and are rejected.
Generated-output filename sanitization and results-directory containment remain
unchanged.

Launch methods wait for the visualizer to close by default, on every platform.
A nonzero process exit is an `IOResult` failure with the exit code and data
filename retained, not a successful launch. Missing executables and Windows
process-start/wait/exit-query failures also return explicit failures. On Windows,
the calling thread suppresses critical loader-error dialogs only around process
creation, restores its previous error-mode flags before waiting, and retains the
original process-start error code. This prevents an invalid executable from
blocking on an interactive dialog in headless CI without changing the
process-wide error mode or relying on the parent's settings. The macOS
Qt app-bundle route reports the status of `/usr/bin/open -W`, not an independently
queried application exit code.

On POSIX, `posix_spawn` receives argument and environment storage prepared entirely
in the parent. There is no MML allocation or environment mutation in a post-`fork`
child. Linux child environments retain inherited variables, prepend the executable's
`lib` directory to `LD_LIBRARY_PATH`, and select its `plugins` and
`plugins/platforms` directories for Qt; the parent's environment is unchanged.
Do not mutate the process environment concurrently with launch, or install a
`SIGCHLD` policy/handler that discards or independently reaps MML's owned children.
Lost child ownership is reported as a wait failure, not success.

A finite timeout (including zero) uses a steady-clock deadline on POSIX.
After expiry, MML sends `SIGTERM`, allows up to 100 ms of grace, then sends
`SIGKILL` if necessary and reaps the owned child, retrying interrupted waits.
Thus cleanup and scheduling can extend elapsed time beyond the requested timeout.
Spawn, wait, signal, and reap errors are explicit failures with OS diagnostics;
a cleanup failure is not reported as an ordinary completed timeout.
The same policy applies to macOS `open -W`: only that owned launcher is terminated
and reaped. LaunchServices may leave the GUI running, including an already-running
application instance. MML neither kills unrelated GUI processes nor claims to know
their exit codes. App-bundle success/failure messages identify the launcher, and
timeout messages explicitly state that the GUI may still be running.

Existing-file lookup tries the input path as supplied, then the results-directory fallback.
Every resolved input is made canonical before process execution, including rigid-body inputs,
so a child's executable-directory working directory does not reinterpret relative data paths.
Empty or null-containing input names, directories supplied as input files, and empty multi-file
lists are rejected before launch. Filesystem failures retain their diagnostics in `IOResult`
instead of being reported as an invalid output filename.

```cpp
const auto invalid = Visualizer::VisualizeRealFunctionFromFile("");
if (!invalid)
    std::cout << invalid.message << '\n'; // No GUI is started.
```

This validation example is also executed by `MML_DocsApp --visualizer-configuration`.
The refactor does not change bundle/backend precedence, default blocking behavior, partial
backend installation support, serialized formats, or vector persistence APIs.

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
    "sine_curve.mml"
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
Visualizer::VisualizeRealFunction(f, "Sample Points", xValues, "samples.mml");
```

---

### Real Functions - Multiple Functions

Visualize multiple functions on the same plot for comparison.

All function and curve collection workflows require nonempty collections with
counts representable as `int`. Pointer collections must contain no null entries.
When legends are supplied, their size must match the collection **exactly**;
neither missing nor extra labels are accepted. This also applies to separate-file
function plots and the per-curve-legend 2D workflow.
The entire collection is checked before creating any output file or launching a
visualizer, so a bad later entry cannot leave earlier output files behind.
Rejection returns `IOError::INVALID_PARAMETERS` with a diagnostic and filename.
Existing file suffixes, titles, sampling, and launch behavior remain unchanged
for valid collections.

Run `MML_DocsApp --collection-contracts` for the bounded, no-GUI validation example
in [`Docs_Demo_CollectionContracts`](../../src/docs_demos/tools/docs_demo_visualizers.cpp).

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
    "trig_comparison.mml"
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
    "interp_comparison.mml"
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
    std::string fileName  // Base name; "_1.mml", "_2.mml" appended
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
    "wave_surface.mml"
);
```

**Output Format:** Grid of (x, y, f(x,y)) values for surface plotting.

---

### Vector Fields

Visualize vector fields F: ℝⁿ → ℝⁿ as arrow plots.

Ordinary 2D and 3D calls export **every sampled vector**, without an implicit
magnitude cutoff. This intentionally removes the former hidden 3D cutoff of
`3`, which could discard an entire field merely because of its units.

Both Cartesian methods also offer an overload with a final `Real upperThreshold`
argument after the filename. Filtering is opt-in: only samples whose L2 magnitude
is **strictly less** than the threshold are exported; samples at the threshold
are omitted, not clamped or rescaled. The threshold is in field-value units and
must be greater than zero; NaN is rejected with `INVALID_PARAMETERS` before
file creation or function evaluation. Positive infinity remains a supported
upper bound, but the ordinary overload is the way to disable filtering entirely.
Pass `3` explicitly to request the former 3D filtering policy.

The runnable `MML_DocsApp --vector-filtering` demo checks both dimensions,
the constant `(4,0,0)` witness, explicit cutoffs, and invalid input without opening
a GUI. Its isolated missing-executable errors are checked separately from the
exported sample counts.

#### 2D Vector Fields

```cpp
static IOResult VisualizeVectorField2DCartesian(
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
    "radial_field.mml"
);
```

**Visualization:** Arrows at each grid point showing direction and magnitude.

#### 3D Vector Fields

```cpp
static IOResult VisualizeVectorField3DCartesian(
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
    "magnetic_field.mml"
);
```

---

**Explicit filtering example:**

```cpp
VectorFunction<3> field([](const VectorN<Real, 3>&) {
    return VectorN<Real, 3>{4, 0, 0};
});
const Real upperThreshold = 3;
auto filtered = Visualizer::VisualizeVectorField3DCartesian(
    field, "Explicit cutoff",
    0, 1, 2, 0, 1, 2, 0, 1, 2,
    "vector_filtered.mml", upperThreshold);
if (!filtered)
    std::cerr << filtered.message << '\n';
```

This explicitly produces zero vector rows because all eight magnitudes are `4`.
Omitting the final argument exports all eight rows. Exporting is followed by
launching the selected visualizer; `IOResult` reports launch failures too.

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
    "circle.mml"
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
    "helix.mml"
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
    "position_time.mml"
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
    "all_states.mml"
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
    "phase_portrait.mml"
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
    "trajectory_3d.mml"
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
Visualizer::VisualizeParticleSimulation2D("collision_sim_100_balls.mml");
```

**Data Format:** File should contain particle positions/velocities over time (format depends on specific visualizer).

#### 3D Particle Simulation

```cpp
static void VisualizeParticleSimulation3D(std::string fileName);
```

**Example:**
```cpp
// 3D N-body simulation
Visualizer::VisualizeParticleSimulation3D("nbody_galaxy.mml");
```

---

## Usage Examples

### Complete Workflow Examples

#### Example 1: Function Analysis

```cpp
#include <mml/tools/Visualizer.h>
#include <mml/base/RealFunction.h>

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
        "gaussian.mml"
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
        "taylor_comparison.mml"
    );
}
```

#### Example 2: Vector Field Visualization

```cpp
#include <mml/tools/Visualizer.h>
#include <mml/base/VectorFunction.h>

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
        "electric_field.mml"
    );
}
```

#### Example 3: ODE Solution Phase Portrait

```cpp
#include <mml/tools/Visualizer.h>
#include <mml/base/ODESystem.h>
#include <mml/algorithms/ODESystemSolver.h>
#include <mml/algorithms/ODEStepCalculators.h>

void visualizePendulum() {
    // Pendulum: θ'' + sin(θ) = 0
    // State: [θ, ω] where ω = θ'
    auto pendulumFunc = [](Real t, const Vector<Real>& state, Vector<Real>& dstate) {
        Real theta = state[0];
        Real omega = state[1];
        dstate[0] = omega;
        dstate[1] = -std::sin(theta);
    };
    
    ODESystem pendulum(2, pendulumFunc);
    Vector<Real> initState({Constants::PI/4.0, 0.0});  // 45° initial angle
    
    // Solve using RK4 fixed step
    ODESystemFixedStepSolver solver(pendulum, StepCalculators::RK4_Basic);
    ODESystemSolution sol = solver.integrate(initState, 0.0, 10.0, 1000);
    
    // Time series of angle
    Visualizer::VisualizeODESysSolCompAsFunc(
        sol, 0, "Angle vs Time", "theta_t.mml"
    );
    
    // Phase portrait (ω vs θ)
    Visualizer::VisualizeODESysSolAsParamCurve2(
        sol, 0, 1, "Phase Portrait", "phase.mml"
    );
    
    // Both variables
    std::vector<std::string> labels = {"θ", "ω"};
    Visualizer::VisualizeODESysSolAsMultiFunc(
        sol, "Pendulum State", labels, "state.mml"
    );
}
```

#### Example 4: 3D Curves

```cpp
#include <mml/tools/Visualizer.h>
#include <mml/core/Curves.h>

void visualize3DCurves() {
    // Trefoil knot (scaled by 50 for visualization)
    auto trefoil = [](Real t) {
        Real x = 50.0 * (std::sin(t) + 2.0*std::sin(2.0*t));
        Real y = 50.0 * (std::cos(t) - 2.0*std::cos(2.0*t));
        Real z = 50.0 * (-std::sin(3.0*t));
        return VectorN<Real, 3>{x, y, z};
    };
    
    ParametricCurve<3> knot(trefoil);
    
    Visualizer::VisualizeParamCurve3D(
        knot,
        "Trefoil Knot",
        0.0, 2.0*Constants::PI,
        500,
        "trefoil.mml"
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

Source code: [MML_Linux_Visualizers](https://github.com/zvanjak/MML_Linux_Visualizers)

#### Qt 6.9 Visualizers
- **RealFuncVisualizer**: 2D function plots
- **SurfaceVisualizer**: 3D surface plots with rotation and lighting
- **ParametricCurve2DVisualizer**: 2D parametric curves
- **ParametricCurve3DVisualizer**: 3D curves with camera controls
- **VectorField2DVisualizer**: 2D arrow plots
- **VectorField3DVisualizer**: 3D vector field visualization
- **Particle2DVisualizer**: 2D particle animations
- **Particle3DVisualizer**: 3D particle/N-body simulations

**Requirements:** Qt 6.9+ runtime libraries.

#### FLTK Visualizers (Lightweight Alternative)
- **RealFuncVisualizer_FLTK**: 2D function plots
- **ParametricCurve2DVisualizer_FLTK**: 2D parametric curves
- **VectorField2DVisualizer_FLTK**: 2D vector field plots

**Requirements:** FLTK 1.3+ (`libfltk1.3-dev` on Debian/Ubuntu).

**When to use FLTK:**
- Minimal dependencies preferred
- Only 2D visualizations needed
- Older systems without Qt 6.9 support

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
- WPF visualizers launched (no external dependencies)
- Paths resolved from the configured bundle root

### Linux

- All visualization methods fully supported
- Qt 6.9 visualizers handle all content (2D and 3D)
- FLTK visualizers available as lightweight alternative for 2D-only workflows
- Requires Qt 6.9+ or FLTK 1.3+ depending on visualizer choice

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
Visualizer::VisualizeRealFunction(f, "Undersampled", 0, 10, 10, "bad.mml");

// Good - captures features
Visualizer::VisualizeRealFunction(f, "Well Sampled", 0, 10, 200, "good.mml");

// Overkill - wastes resources
Visualizer::VisualizeRealFunction(f, "Oversampled", 0, 10, 10000, "slow.mml");
```

**Guidelines:**
- Smooth functions: 100-200 points
- Oscillatory functions: 10-20 points per period
- Discontinuities: Dense sampling near jumps

### 2. Use Descriptive Titles and Filenames

```cpp
// Poor
Visualizer::VisualizeRealFunction(f, "f", 0, 1, 100, "out.mml");

// Good
Visualizer::VisualizeRealFunction(
    f, 
    "Velocity Profile (Re=1000)", 
    0, 1, 100, 
    "velocity_profile_re1000.mml"
);
```

### 3. Separate Files for Different Visualizations

```cpp
// Organize results by category
Visualizer::VisualizeRealFunction(f, "...", 0, 1, 100, "analysis/case1_func.mml");
Visualizer::VisualizeScalarFunc2DCartesian(s, "...", -1, 1, 50, -1, 1, 50, "analysis/case1_field.mml");
```

### 4. Verify Data Before Visualization

```cpp
// Check ODE solution validity before visualizing
if (sol.IsValid()) {
    Visualizer::VisualizeODESysSolAsMultiFunc(sol, "...", legend, "ode_sol.mml");
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
    std::string filename = "method_" + methods[i] + ".mml";
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
Visualizer::VisualizeRealFunction(f, "...", 0, 1, 100, "data.mml");
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
    surface, "Dense Grid", -10, 10, nx, -10, 10, ny, "dense.mml"
);
// Requires ~30 MB file + evaluation overhead
```

**Recommendation:** Start with coarse sampling, refine as needed.

---

## Integration with Other Components

### With Serializer

```cpp
// Manual export for custom processing
Serializer::SaveRealFunc(f, "Function", 0, 10, 100, "results/data.mml");

// Then use Visualizer for same data (re-exports)
Visualizer::VisualizeRealFunction(f, "Function", 0, 10, 100, "data.mml");
```

**Tip:** Use `Serializer` directly when you only need data export without visualization.

### With ODE Solvers

```cpp
// Solve system
ODESystemFixedStepSolver<4> solver(system);
auto sol = solver.integrate_RK4(initState, 0.0, 100.0, 0.01);

// Visualize multiple aspects
Visualizer::VisualizeODESysSolAsMultiFunc(sol, "All Vars", labels, "all.mml");
Visualizer::VisualizeODESysSolAsParamCurve2(sol, 0, 1, "Phase", "phase.mml");
Visualizer::VisualizeODESysSolCompAsFunc(sol, 0, "X vs t", "x_t.mml");
```

### With Function Analysis

```cpp
#include <mml/algorithms/Analyzers/FunctionsAnalyzer.h>

// Find and visualize extrema
auto criticalPoints = FunctionsAnalyzer::FindCriticalPoints(f, a, b, dx);

// Visualize function with critical points marked
Vector<Real> xCritical(criticalPoints.size());
for (size_t i = 0; i < criticalPoints.size(); ++i) {
    xCritical[i] = criticalPoints[i].position;
}

Visualizer::VisualizeRealFunction(f, "Function", a, b, 200, "func.mml");
Visualizer::VisualizeRealFunction(f, "Critical Points", xCritical, "critical.mml");
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
// Check bundle configuration in MMLVisualizators.h
// Verify GetRealFuncVisualizerPath() returns valid path

// Data is still exported - use external tool:
// Python: plt.loadtxt("results/data.mml"); plt.plot(...)
// gnuplot: plot "results/data.mml" using 1:2 with lines
// MATLAB: data = load("results/data.mml"); plot(data(:,1), data(:,2))
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
Visualizer::VisualizeRealFunction(f, "Better", 0, 10, 500, "better.mml");  // Was 50

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
