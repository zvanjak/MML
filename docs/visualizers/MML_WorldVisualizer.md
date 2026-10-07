# MML_WorldVisualizer

A general-purpose 3D scene visualizer for MML. While the other MML visualizers each display a single type of data (a function graph, a vector field, particle trajectories, etc.), `MML_WorldVisualizer` lets you compose arbitrary 3D scenes by placing multiple primitives — points, vectors, lines, planes, patches — into a single view. It is the natural tool for geometric illustrations, differential-forms visualizations, field glyphs, and any scene that mixes object types.

## How It Works

You write a plain-text `.mmlworld` scene file, then launch `MML_WorldVisualizer` with that file as an argument. The visualizer parses the file line by line, constructs 3D meshes for each primitive, and renders them in an interactive WPF viewport with mouse-driven camera control.

```
[MML C++ code] → write .mmlworld file → Visualizer::VisualizeWorldSceneFromFile() → MML_WorldVisualizer.exe
```

Scene files are saved by convention to the results folder (see `GetResultFilesPath()`), and the `Visualizer` class resolves bare filenames against that folder automatically.

---

## Launching from MML C++ Code

```cpp
#include <mml/tools/Visualizer.h>

// Write your scene to a file, then launch:
VisualizerResult result = Visualizer::VisualizeWorldSceneFromFile("my_scene.mmlworld");

if (!result.success)
    std::cerr << "Visualizer error: " << result.errorMessage << "\n";
```

You can also pass the full path of a file written anywhere:

```cpp
VisualizerResult result = Visualizer::VisualizeWorldSceneFromFile("/tmp/debug_scene.mmlworld");
```

**Environment variable**: set `MML_VISUALIZER_SAVE_ONLY=1` to skip launching the GUI (useful in CI or headless environments). The scene file is still written.

---

## Scene File Format

Scene files are UTF-8 plain text. Each non-empty, non-comment line is one command. Lines starting with `#` are comments and are ignored.

### Header

```
MML_WORLD_SCENE 1
```

The first non-comment line must be exactly this. The `1` is the format version.

### Optional setup commands

| Command | Syntax | Description |
|---------|--------|-------------|
| `TITLE` | `TITLE <text>` | Sets the window title bar text |
| `CAMERA` | `CAMERA <x> <y> <z>` | Sets initial camera position |

---

## Supported Primitives

### `POINT`
Renders a sphere at a given location.

```
POINT <px> <py> <pz>  <radius>  <color>
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `px py pz` | float | Center position |
| `radius` | float | Sphere radius |
| `color` | hex | `#RRGGBB` |

**Example:**
```
POINT 0 0 0  3.4  #FFFFFF
POINT 50 50 50  8  #FF4444
```

---

### `VECTOR`
Renders an arrow (cylinder body + cone tip) representing a vector rooted at a position.

```
VECTOR <px> <py> <pz>  <vx> <vy> <vz>  <radius>  <color>
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `px py pz` | float | Arrow base (origin of vector) |
| `vx vy vz` | float | Vector components (determines direction **and** length) |
| `radius` | float | Shaft radius |
| `color` | hex | `#RRGGBB` |

The arrow length is `|v| * scale` — you control the visual scale by pre-multiplying the vector before writing. A zero-length vector is silently skipped.

**Example:**
```
VECTOR 0 0 0  100 100 100  4.0  #FFD23F
VECTOR 0 0 0   50   0   0  2.5  #4CB3FF
```

---

### `LINE`
Renders a thin cylinder between two 3D points.

```
LINE <x0> <y0> <z0>  <x1> <y1> <z1>  <radius>  <color>
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `x0 y0 z0` | float | Start point |
| `x1 y1 z1` | float | End point |
| `radius` | float | Tube radius |
| `color` | hex | `#RRGGBB` |

Useful for coordinate axes, field-line segments, grid lines, and polygon edges.

**Example:**
```
LINE -160 0 0  160 0 0  1.2  #777777
LINE 0 -160 0  0 160 0  1.2  #777777
LINE 0 0 -80   0 0 180  1.2  #777777
```

---

### `PLANE`
Renders a square plane perpendicular to a given normal vector, centered at a point.

```
PLANE <cx> <cy> <cz>  <nx> <ny> <nz>  <size>  <color>  <opacity>
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `cx cy cz` | float | Center of the plane |
| `nx ny nz` | float | Normal vector (direction the plane faces) |
| `size` | float | Side length of the square |
| `color` | hex | `#RRGGBB` |
| `opacity` | float | 0.0 (transparent) … 1.0 (opaque) |

Rendered double-sided. A zero-length normal is silently skipped.

**Example — stacked planes (one-form glyph):**
```
PLANE -12 -6 -18  0.535 -0.267 0.802  145.0  #56D6C9  0.22
PLANE   0  0   0  0.535 -0.267 0.802  145.0  #56D6C9  0.22
PLANE  12  6  18  0.535 -0.267 0.802  145.0  #56D6C9  0.22
```

---

### `PATCH`
Renders a parallelogram (quad) defined by two edge vectors from a center point.

```
PATCH <cx> <cy> <cz>  <ux> <uy> <uz>  <vx> <vy> <vz>  <color>  <opacity>
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `cx cy cz` | float | Center of the patch |
| `ux uy uz` | float | First edge vector (half-width in u direction) |
| `vx vy vz` | float | Second edge vector (half-width in v direction) |
| `color` | hex | `#RRGGBB` |
| `opacity` | float | 0.0 … 1.0 |

The four corners are: `center ± 0.5u ± 0.5v`. Rendered double-sided.

**Example — flux two-form glyph:**
```
PATCH 0 0 0  80 0 0  0 60 0  #FF6B35  0.50
```

---

## Color Specification

All color parameters use hexadecimal color strings:

| Format | Example | Description |
|--------|---------|-------------|
| `#RRGGBB` | `#FF6B35` | Opaque color (alpha = 255) |
| `#AARRGGBB` | `#80FF6B35` | Color with explicit alpha channel |

When `opacity` is given as a separate parameter (for `PLANE` and `PATCH`), it multiplies the brush opacity on top of any alpha embedded in the color.

---

## Comments and Blank Lines

```
# This is a comment - ignored by the parser
# Use comments to document your scene format inline

MML_WORLD_SCENE 1
TITLE My Scene

# VECTOR px py pz  vx vy vz  radius color
VECTOR 0 0 0  100 0 0  3.0  #FF0000
```

---

## Complete Example

The following generates a scene with a coordinate frame, a vector, its dual one-form planes, and a flux patch — the canonical "typed forms" visualization:

```cpp
#include <mml/tools/Visualizer.h>
#include <fstream>
#include <iomanip>
#include <filesystem>

void WriteDemoScene()
{
    // Scene file goes to the standard results folder
    std::filesystem::path results = GetResultFilesPath();
    std::filesystem::create_directories(results);
    std::string path = (results / "demo.mmlworld").string();

    std::ofstream scene(path);
    scene << std::fixed << std::setprecision(4);

    scene << "MML_WORLD_SCENE 1\n";
    scene << "TITLE Demo: vector + one-form planes + flux patch\n";
    scene << "CAMERA 420 260 280\n\n";

    // Coordinate axes
    scene << "LINE -150 0 0  150 0 0  1.2  #666666\n";
    scene << "LINE 0 -150 0  0 150 0  1.2  #666666\n";
    scene << "LINE 0 0 -80   0 0 180  1.2  #666666\n\n";

    // Vector v = (2, -1, 3), scaled by 35
    scene << "VECTOR 0 0 0   70 -35 105   4.0  #FFD23F\n\n";

    // One-form planes perpendicular to v (stacked at -1, 0, +1 levels)
    double nx = 0.5345, ny = -0.2673, nz = 0.8018; // unit(v)
    for (int i = -1; i <= 1; i++) {
        double cx = nx * 22.0 * i;
        double cy = ny * 22.0 * i;
        double cz = nz * 22.0 * i;
        scene << "PLANE " << cx << " " << cy << " " << cz
              << "  " << nx << " " << ny << " " << nz
              << "  145.0  #56D6C9  0.22\n";
    }
    scene << "\n";

    // Flux two-form patch
    scene << "PATCH 0 0 0   80 0 0   0 60 0   #FF6B35  0.50\n\n";

    // Origin marker
    scene << "POINT 0 0 0  3.5  #FFFFFF\n";

    scene.close();
    Visualizer::VisualizeWorldSceneFromFile(path);
}
```

---

## Interactive Controls

| Action | Effect |
|--------|--------|
| Left-drag | Rotate scene (trackball) |
| Right-drag | Pan camera |
| Scroll wheel | Zoom in / out |
| **Look at Center** button | Reset view to origin |
| **Load Data** button | Open a `.mmlworld` file via file dialog |
| **Export Image** button | Save current viewport as PNG/JPEG/BMP |
| **Toggle Theme** button | Switch between dark and light background |

---

## Typical Use Cases

| Goal | Primitives Used |
|------|----------------|
| Visualize a vector field at sample points | `VECTOR` + `POINT` |
| Illustrate a differential one-form | `VECTOR` + stacked `PLANE` |
| Illustrate a differential two-form | `PATCH` + `LINE` |
| Combined typed-forms glyph | `VECTOR` + `PLANE` + `PATCH` + `POINT` |
| Show coordinate frame / basis | `LINE` + `VECTOR` |
| Mark specific points in a scene | `POINT` |
| Draw polygon edges or segments | `LINE` |

---

## Related Files

| File | Description |
|------|-------------|
| [mml/tools/Visualizer.h](../../mml/tools/Visualizer.h) | C++ API — `Visualizer::VisualizeWorldSceneFromFile()` |
| [mml/MMLVisualizators.h](../../mml/MMLVisualizators.h) | Path helper — `GetWorldVisualizerPath()` |
| [src/visualization_examples/show_typed_forms_3d.cpp](../../src/visualization_examples/show_typed_forms_3d.cpp) | Full working example (vortex tube + magnetic dipole field scenes) |
| [tools/visualizers/win/WPF/MML_WorldVisualizer/](../../tools/visualizers/win/WPF/MML_WorldVisualizer/) | Bundled pre-built visualizer app |
| `D:\Projects\MML_Visualizers\WPF\MML_WorldVisualizer\` | Source code of the WPF application |
