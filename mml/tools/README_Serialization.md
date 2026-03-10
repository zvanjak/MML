# MML Serialization & Data I/O Guide

> **Comprehensive documentation for MinimalMathLibrary's serialization, data loading, and file I/O capabilities.**

---

## Table of Contents

1. [Overview](#overview)
2. [Quick Reference](#quick-reference)
3. [Matrix I/O](#matrix-io)
4. [Vector I/O](#vector-io)
5. [DataLoader (CSV, TSV, JSON)](#dataloader)
6. [Serializer (Visualization Output)](#serializer)
7. [RigidBody Simulation Serializer](#rigidbody-simulation-serializer)
8. [File Format Specifications](#file-format-specifications)
9. [Best Practices](#best-practices)
10. [Examples](#examples)

---

## Overview

MML provides three main serialization subsystems:

| Subsystem | Purpose | Direction | Header |
|-----------|---------|-----------|--------|
| **Matrix/Vector I/O** | Save/load basic linear algebra objects | Bidirectional | `<MML/base/Matrix/Matrix.h>`, `<MML/base/Vector/Vector.h>` |
| **DataLoader** | Load tabular data from files | Read-only | `<MML/core/DataLoader.h>` |
| **Serializer** | Export mathematical objects for visualization | Write-only | `<MML/core/Serializer.h>` |

### Supported File Formats

| Format | Extension | Binary? | Description |
|--------|-----------|---------|-------------|
| MML Text | `.mml` | No | Human-readable visualization format |
| MML Matrix Binary | `.mmlm` | Yes | Efficient matrix storage |
| MML Vector Binary | `.mmlv` | Yes | Efficient vector storage |
| CSV | `.csv` | No | Comma-separated values |
| TSV | `.tsv`, `.tab` | No | Tab-separated values |
| JSON | `.json` | No | Array of flat objects |

---

## Quick Reference

### Matrix Operations
```cpp
#include <MML/base/Matrix/Matrix.h>

Matrix<Real> A(100, 100);

// Text format (human-readable)
A.SaveToFile("matrix.txt");
Matrix<Real>::LoadFromFile("matrix.txt", A);

// CSV format
A.SaveToCSV("matrix.csv");
Matrix<Real>::LoadFromCSV("matrix.csv", A);

// Binary format (fastest, smallest)
A.SaveToBinaryFile("matrix.mmlm");
Matrix<Real>::LoadFromBinaryFile("matrix.mmlm", A);
```

### Vector Operations
```cpp
#include <MML/base/Vector/Vector.h>

Vector<Real> v(1000);

// Text format
v.SaveToFile("vector.txt");
Vector<Real>::LoadFromFile("vector.txt", v);

// Binary format
v.SaveToBinaryFile("vector.mmlv");
Vector<Real>::LoadFromBinaryFile("vector.mmlv", v);
```

### DataLoader Operations
```cpp
#include <MML/core/DataLoader.h>
using namespace MML::Data;

// Auto-detect format by extension
auto dataset = LoadData("data.csv");  // or .tsv, .json

// Explicit format
auto csv_data = LoadCSV("data.csv");
auto json_data = LoadJSON("data.json");

// Access data
for (const auto& row : dataset.rows) {
    Real x = row.GetReal("x");
    std::string name = row.GetString("name");
}
```

### Serializer Operations
```cpp
#include <MML/core/Serializer.h>
using namespace MML::Serializer;

// Save function evaluations
RealFunction f = [](Real x) { return std::sin(x); };
SaveRealFunction(f, "sin.mml", 0.0, 6.28, 100);

// Save ODE solution
SaveODESolution(solution, "ode.mml");

// Save parametric curve
SaveParametricCurve3D(curve, "helix.mml", 0.0, 10.0, 500);
```

---

## Matrix I/O

### Header
```cpp
#include <MML/base/Matrix/Matrix.h>
```

### Save Methods

#### `SaveToFile(filename)` - Text Format
Saves matrix in space-separated text format with dimensions header.

```cpp
Matrix<Real> A(3, 4);
// ... fill matrix ...
A.SaveToFile("matrix.txt");
```

**Output format:**
```
3 4
1.234567890123456 2.345678901234567 3.456789012345678 4.567890123456789
5.678901234567890 6.789012345678901 7.890123456789012 8.901234567890123
9.012345678901234 0.123456789012345 1.234567890123456 2.345678901234567
```

#### `SaveToCSV(filename, precision)` - CSV Format
Saves as comma-separated values with optional precision control.

```cpp
A.SaveToCSV("matrix.csv");           // Default precision (15)
A.SaveToCSV("matrix.csv", 6);        // 6 decimal places
```

**Output format:**
```csv
1.234568,2.345679,3.456789,4.567890
5.678901,6.789012,7.890123,8.901234
9.012346,0.123457,1.234568,2.345679
```

#### `SaveToBinaryFile(filename)` - Binary Format
Most efficient storage. Uses magic number `MMLM` for format identification.

```cpp
A.SaveToBinaryFile("matrix.mmlm");
```

**Binary structure:**
```
[4 bytes] Magic: "MMLM"
[4 bytes] Rows (uint32)
[4 bytes] Cols (uint32)
[R×C × 8 bytes] Data (row-major, double)
```

### Load Methods

#### `LoadFromFile(filename, matrix)` - Text Format
```cpp
Matrix<Real> A;
Matrix<Real>::LoadFromFile("matrix.txt", A);
// A is now resized and filled
```

#### `LoadFromCSV(filename, matrix)` - CSV Format
```cpp
Matrix<Real> A;
Matrix<Real>::LoadFromCSV("matrix.csv", A);
```

#### `LoadFromBinaryFile(filename, matrix)` - Binary Format
```cpp
Matrix<Real> A;
Matrix<Real>::LoadFromBinaryFile("matrix.mmlm", A);
```

### Performance Comparison

| Format | 1000×1000 Matrix | File Size | Read Speed |
|--------|------------------|-----------|------------|
| Text | ~23 MB | Largest | Slowest |
| CSV | ~18 MB | Medium | Medium |
| Binary | ~8 MB | Smallest | Fastest |

**Recommendation:** Use binary for large matrices, text/CSV for debugging and interoperability.

---

## Vector I/O

### Header
```cpp
#include <MML/base/Vector/Vector.h>
```

### Save Methods

#### `SaveToFile(filename)` - Text Format
Saves one value per line with dimension header.

```cpp
Vector<Real> v(5);
v[0] = 1.1; v[1] = 2.2; v[2] = 3.3; v[3] = 4.4; v[4] = 5.5;
v.SaveToFile("vector.txt");
```

**Output format:**
```
5
1.100000000000000
2.200000000000000
3.300000000000000
4.400000000000000
5.500000000000000
```

#### `SaveToBinaryFile(filename)` - Binary Format
Uses magic number `MMLV`.

```cpp
v.SaveToBinaryFile("vector.mmlv");
```

**Binary structure:**
```
[4 bytes] Magic: "MMLV"
[4 bytes] Size (uint32)
[N × 8 bytes] Data (double)
```

### Load Methods

```cpp
Vector<Real> v;
Vector<Real>::LoadFromFile("vector.txt", v);
Vector<Real>::LoadFromBinaryFile("vector.mmlv", v);
```

---

## DataLoader

### Header
```cpp
#include <MML/core/DataLoader.h>
using namespace MML::Data;
```

### Core Types

```cpp
// A single data value (type-safe variant)
struct DataValue {
    enum class Type { Null, Bool, Integer, Real, String };
    
    bool IsNull() const;
    bool GetBool() const;
    int64_t GetInt() const;
    double GetReal() const;
    std::string GetString() const;
};

// A row of named values
struct DataRow {
    DataValue Get(const std::string& column) const;
    Real GetReal(const std::string& column) const;
    std::string GetString(const std::string& column) const;
    // ... convenience methods
};

// Complete dataset
struct Dataset {
    std::vector<std::string> columns;
    std::vector<DataRow> rows;
    
    size_t NumRows() const;
    size_t NumColumns() const;
    bool HasColumn(const std::string& name) const;
};
```

### Loading Functions

#### Auto-Detect by Extension
```cpp
Dataset data = LoadData("file.csv");   // Uses LoadCSV
Dataset data = LoadData("file.tsv");   // Uses LoadTSV  
Dataset data = LoadData("file.json");  // Uses LoadJSON
```

#### Explicit Loaders
```cpp
Dataset csv_data = LoadCSV("data.csv");
Dataset tsv_data = LoadTSV("data.tsv");
Dataset json_data = LoadJSON("data.json");

// With options
LoadOptions opts;
opts.hasHeader = true;
opts.delimiter = ';';  // For European CSV
Dataset data = LoadCSV("european.csv", opts);
```

### CSV/TSV Format Requirements

**With header (default):**
```csv
x,y,z,label
1.0,2.0,3.0,point1
4.0,5.0,6.0,point2
```

**Without header:**
```cpp
LoadOptions opts;
opts.hasHeader = false;
opts.columnNames = {"x", "y", "z", "label"};
auto data = LoadCSV("noheader.csv", opts);
```

### JSON Format Requirements

The JSON parser expects an **array of flat objects**:

```json
[
    {"x": 1.0, "y": 2.0, "label": "A"},
    {"x": 3.0, "y": 4.0, "label": "B"},
    {"x": 5.0, "y": 6.0, "label": "C"}
]
```

**Supported value types:**
- Numbers (integers and floating-point)
- Strings (with escape sequences)
- Booleans (`true`, `false`)
- Null (`null`)

**Not supported:**
- Nested objects
- Nested arrays
- Top-level objects (must be array)

### Accessing Data

```cpp
Dataset data = LoadCSV("points.csv");

// Iterate all rows
for (const auto& row : data.rows) {
    Real x = row.GetReal("x");
    Real y = row.GetReal("y");
    std::string label = row.GetString("label");
    
    // Or use type-safe access
    DataValue val = row.Get("x");
    if (!val.IsNull()) {
        double x = val.GetReal();
    }
}

// Check columns
if (data.HasColumn("weight")) {
    // process weight column
}

// Get dimensions
std::cout << "Loaded " << data.NumRows() << " rows, " 
          << data.NumColumns() << " columns\n";
```

---

## Serializer

The Serializer namespace provides **write-only** output for visualization tools.

### Header
```cpp
#include <MML/core/Serializer.h>
using namespace MML::Serializer;
```

### Supported Object Types

| Object Type | Function | Output Format |
|-------------|----------|---------------|
| Real functions | `SaveRealFunction()` | `REAL_FUNCTION_SAMPLED` |
| Multiple functions | `SaveMultiFunctions()` | `MULTI_REAL_FUNCTION_SAMPLED` |
| Parametric curve 2D | `SaveParametricCurve2D()` | `PARAMETRIC_CURVE_CARTESIAN_2D` |
| Parametric curve 3D | `SaveParametricCurve3D()` | `PARAMETRIC_CURVE_CARTESIAN_3D` |
| Parametric surface | `SaveParametricSurface()` | `PARAMETRIC_SURFACE_CARTESIAN_3D` |
| Scalar function 2D | `SaveScalarFunction2D()` | `SCALAR_FUNCTION_CARTESIAN_2D` |
| Scalar function 3D | `SaveScalarFunction3D()` | `SCALAR_FUNCTION_CARTESIAN_3D` |
| Vector field 2D | `SaveVectorField2D()` | `VECTOR_FIELD_CARTESIAN_2D` |
| Vector field 3D | `SaveVectorField3D()` | `VECTOR_FIELD_CARTESIAN_3D` |
| ODE solution (single) | `SaveODESolution()` | `MULTI_REAL_FUNCTION_VARIABLE_SAMPLED` |
| ODE solution (system) | `SaveODESystemSolution()` | `MULTI_REAL_FUNCTION_VARIABLE_SAMPLED` |
| Particle simulation 2D | `SaveParticleSimulation2D()` | `PARTICLE_SIMULATION_2D` |
| Particle simulation 3D | `SaveParticleSimulation3D()` | `PARTICLE_SIMULATION_3D` |
| Field lines 2D | `SaveFieldLines2D()` | `FIELD_LINES_2D` |
| Field lines 3D | `SaveFieldLines3D()` | `FIELD_LINES_3D` |

### Examples

#### Save a Real Function
```cpp
auto f = [](Real x) { return std::sin(x) * std::exp(-x/10); };

SaveRealFunction(f, "damped_sine.mml", 
    0.0,    // x_start
    20.0,   // x_end  
    200     // num_points
);
```

#### Save Multiple Functions Together
```cpp
std::vector<std::function<Real(Real)>> funcs = {
    [](Real x) { return std::sin(x); },
    [](Real x) { return std::cos(x); },
    [](Real x) { return std::sin(2*x); }
};
std::vector<std::string> names = {"sin", "cos", "sin2x"};

SaveMultiFunctions(funcs, names, "trig.mml", 0.0, 6.28, 100);
```

#### Save Parametric Curve 3D (Helix)
```cpp
auto helix = [](Real t) {
    return Vector3Cartesian(
        std::cos(t),      // x
        std::sin(t),      // y
        t / (2 * M_PI)    // z
    );
};

SaveParametricCurve3D(helix, "helix.mml", 
    0.0,     // t_start
    20.0,    // t_end
    500      // num_points
);
```

#### Save ODE Solution
```cpp
ODESystem system = ...;
ODESolver solver;
auto solution = solver.Solve(system, y0, t_span);

SaveODESolution(solution, "pendulum.mml");
```

#### Save Scalar Function 2D (Heat Map)
```cpp
auto f = [](Real x, Real y) { return std::sin(x) * std::cos(y); };

SaveScalarFunction2D(f, "sincos.mml",
    -3.14, 3.14,  // x range
    -3.14, 3.14,  // y range
    50, 50        // grid resolution
);
```

#### Save Vector Field 3D
```cpp
auto field = [](const Vector3Cartesian& p) {
    return Vector3Cartesian(-p.Y(), p.X(), 0);  // Rotation field
};

SaveVectorField3D(field, "rotation.mml",
    -2, 2, -2, 2, -2, 2,  // bounds
    10, 10, 10            // grid resolution
);
```

---

## RigidBody Simulation Serializer

Specialized serializer for rigid body physics simulations.

### Header
```cpp
#include "RigidBody/RigidBodySerializer.h"
using namespace MPL;
```

### Two Output Modes

#### Mode 1: Single File (Recommended)
All bodies in one file with `RIGID_BODY_SIMULATION` format.

```cpp
RigidBodySerializer serializer;
serializer.SetContainerSize(5.0);  // 10m cube container

// Add body metadata
serializer.AddBody(RigidBodyInfo::Box(10.0, 1.0, 0.5, 0.3, "Red", "Box1"));
serializer.AddBody(RigidBodyInfo::Box(8.0, 0.75, 0.5, 0.4, "Blue", "Box2"));
serializer.AddBody(RigidBodyInfo::Sphere(5.0, 0.5, "LimeGreen", "Sphere1"));

// Run simulation...
simulator.Run();

// Save all bodies to one file
auto result = serializer.SaveAllBodies(
    "simulation.mml",
    simulator.GetHistory(),
    config.timeStep,
    10  // Save every 10th frame
);

if (!result) {
    std::cerr << "Error: " << result.message << "\n";
}
```

#### Mode 2: Separate Files (Legacy)
One file per body with `RIGID_BODY_TRAJECTORY_3D` format.

```cpp
serializer.SaveSingleBody("box1.mml", history, 0, dt, 10);
serializer.SaveSingleBody("box2.mml", history, 1, dt, 10);
serializer.SaveSingleBody("sphere.mml", history, 2, dt, 10);
```

### Visualization

Both formats are supported by the WPF visualizer:

```cpp
// Launch visualizer automatically
Visualizer::VisualizeRigidBodySimulation("simulation.mml");

// Or with multiple separate files
std::vector<std::string> files = {"box1.mml", "box2.mml", "sphere.mml"};
Visualizer::VisualizeRigidBodyTrajectory(files);
```

---

## File Format Specifications

### MML Text Header Format

All `.mml` visualization files start with:
```
# MML Visualization Data File
# Format: <FORMAT_TYPE>
# Generated: <timestamp>
# 
# <format-specific metadata>
#
<data section>
```

### REAL_FUNCTION_SAMPLED
```
# Format: REAL_FUNCTION_SAMPLED
# Points: 100
#
# x, f(x)
0.000000000000000 0.000000000000000
0.063466518254339 0.063424925152302
...
```

### PARAMETRIC_CURVE_CARTESIAN_3D
```
# Format: PARAMETRIC_CURVE_CARTESIAN_3D
# Points: 500
#
# t, x, y, z
0.000000 1.000000 0.000000 0.000000
0.040080 0.999196 0.040067 0.006380
...
```

### RIGID_BODY_SIMULATION
```
# MML Rigid Body Simulation
# Format: RIGID_BODY_SIMULATION
# Version: 1.0
#
CONTAINER_SIZE: 5.0

BODIES: 3
BODY 0: BOX mass=10.0 half_extents=(1.0, 0.5, 0.3) color=Red name=Box1
BODY 1: BOX mass=8.0 half_extents=(0.75, 0.5, 0.4) color=Blue name=Box2
BODY 2: SPHERE mass=5.0 radius=0.5 color=LimeGreen name=Sphere1

FRAMES: 401
FRAME 0 t=0.000000
  0: pos=(x,y,z) quat=(w,x,y,z) vel=(vx,vy,vz) omega=(wx,wy,wz)
  1: pos=(x,y,z) quat=(w,x,y,z) vel=(vx,vy,vz) omega=(wx,wy,wz)
  2: pos=(x,y,z) quat=(w,x,y,z) vel=(vx,vy,vz) omega=(wx,wy,wz)
FRAME 1 t=0.010000
...
```

### Binary Format Constants

All magic numbers and version constants are centralized in `MML::BinaryFormat` namespace (`MMLBase.h`):

```cpp
namespace MML::BinaryFormat {
    // Magic numbers (4-byte ASCII identifiers)
    constexpr uint32_t MAGIC_MATRIX = 0x4D4D4C4D;  // "MMLM"
    constexpr uint32_t MAGIC_VECTOR = 0x4D4D4C56;  // "MMLV"
    constexpr uint32_t MAGIC_SPARSE = 0x4D4D4C53;  // "MMLS" (reserved)
    constexpr uint32_t MAGIC_TENSOR = 0x4D4D4C54;  // "MMLT" (reserved)
    
    // Current format versions
    constexpr uint32_t VERSION_MATRIX = 1;
    constexpr uint32_t VERSION_VECTOR = 1;
    
    // File extensions
    constexpr const char* EXT_MATRIX = "mmlm";
    constexpr const char* EXT_VECTOR = "mmlv";
}
```

### Binary Matrix Format (`.mmlm`)
```
Offset  Size    Content
0       4       Magic: BinaryFormat::MAGIC_MATRIX (0x4D4D4C4D = "MMLM")
4       4       Version: BinaryFormat::VERSION_MATRIX (currently 1)
8       4       Rows (uint32, little-endian)
12      4       Columns (uint32, little-endian)
16      4       Element size in bytes (sizeof(Type))
20      4       Reserved (0)
24      R×C×8   Data (double, row-major, little-endian)
```

### Binary Vector Format (`.mmlv`)
```
Offset  Size    Content
0       4       Magic: BinaryFormat::MAGIC_VECTOR (0x4D4D4C56 = "MMLV")
4       4       Version: BinaryFormat::VERSION_VECTOR (currently 1)
8       4       Size (uint32, little-endian)
12      4       Element size in bytes (sizeof(Type))
16      N×8     Data (double, little-endian)
```

---

## Best Practices

### 1. Choose the Right Format

| Scenario | Recommended Format |
|----------|--------------------|
| Large matrices (>1000×1000) | Binary `.mmlm` |
| Debugging / inspection | Text or CSV |
| Interoperability (Excel, Python) | CSV |
| Visualization | `.mml` via Serializer |
| Simulation playback | `.mml` RigidBodySerializer |

### 2. Handle Errors

```cpp
// Matrix I/O returns bool
if (!A.SaveToFile("matrix.txt")) {
    std::cerr << "Failed to save matrix\n";
}

Matrix<Real> B;
if (!Matrix<Real>::LoadFromFile("matrix.txt", B)) {
    std::cerr << "Failed to load matrix\n";
}

// RigidBodySerializer returns SerializeResult
auto result = serializer.SaveAllBodies(...);
if (!result) {
    std::cerr << "Error: " << result.message << "\n";
}
```

### 3. Use Appropriate Precision

```cpp
// For CSV, control decimal places
A.SaveToCSV("matrix.csv", 6);   // 6 decimals (smaller file)
A.SaveToCSV("matrix.csv", 15);  // Full precision (default)

// Binary always uses full double precision
A.SaveToBinaryFile("matrix.mmlm");  // Always 15+ significant digits
```

### 4. Downsample Large Simulations

```cpp
// Save every Nth frame to reduce file size
int skipFactor = 10;  // 1000 Hz → 100 Hz
serializer.SaveAllBodies(filename, history, dt, skipFactor);
```

### 5. Validate Before Processing

```cpp
Dataset data = LoadCSV("input.csv");

// Check expected columns exist
if (!data.HasColumn("x") || !data.HasColumn("y")) {
    throw std::runtime_error("Missing required columns");
}

// Check row count
if (data.NumRows() == 0) {
    throw std::runtime_error("Empty dataset");
}
```

---

## Examples

### Complete Example: Load CSV → Process → Save Binary

```cpp
#include <MML/base/Matrix/Matrix.h>
#include <MML/core/DataLoader.h>
using namespace MML::Data;

int main() {
    // Load data from CSV
    Dataset data = LoadCSV("measurements.csv");
    
    // Convert to matrix manually (extract numeric columns)
    size_t n = data.NumRows();
    Matrix<Real> A(n, 3);
    
    for (size_t i = 0; i < n; ++i) {
        A(i, 0) = data.rows[i].GetReal("x");
        A(i, 1) = data.rows[i].GetReal("y");
        A(i, 2) = data.rows[i].GetReal("z");
    }
    
    // Process...
    Matrix<Real> result = SomeComputation(A);
    
    // Save efficiently
    result.SaveToBinaryFile("result.mmlm");
    
    return 0;
}
```

### Complete Example: Serialize Function for Visualization

```cpp
#include <MML/core/Serializer.h>
#include <MML/core/Visualizer.h>
using namespace MML::Serializer;

int main() {
    // Define a function
    auto bessel = [](Real x) { return std::cyl_bessel_j(0, x); };
    
    // Save for visualization
    SaveRealFunction(bessel, "bessel_j0.mml", 0.0, 20.0, 500);
    
    // Launch visualizer
    Visualizer::VisualizeRealFunction("bessel_j0.mml");
    
    return 0;
}
```

### Complete Example: Rigid Body Simulation Pipeline

```cpp
#include "RigidBody/RigidBodies.h"
#include "RigidBody/RigidBodySimulator.h"
#include "RigidBody/RigidBodySerializer.h"
#include <MML/core/Visualizer.h>

int main() {
    // Create bodies
    RigidBodyBox box(10.0, 1.0, 0.5, 0.3);
    box.Position() = Vec3Cart(0, 0, 2);
    box.Velocity() = Vec3Cart(1, 0, -1);
    
    RigidBodySphere sphere(5.0, 0.5);
    sphere.Position() = Vec3Cart(2, 0, 3);
    sphere.Velocity() = Vec3Cart(-1, 0, -2);
    
    // Configure simulation
    SimulationConfig config;
    config.timeStep = 0.001;
    config.totalTime = 10.0;
    config.containerHalfSize = 5.0;
    config.coeffRestitution = 1.0;
    
    // Run simulation
    RigidBodySimulator sim(config);
    sim.AddBody(box);
    sim.AddBody(sphere);
    sim.Run();
    
    // Serialize
    RigidBodySerializer serializer;
    serializer.SetContainerSize(config.containerHalfSize);
    serializer.AddBody(RigidBodyInfo::Box(10.0, 1.0, 0.5, 0.3, "Red", "Box"));
    serializer.AddBody(RigidBodyInfo::Sphere(5.0, 0.5, "Green", "Sphere"));
    
    serializer.SaveAllBodies("simulation.mml", sim.GetHistory(), 
                             config.timeStep, 10);
    
    // Visualize
    Visualizer::VisualizeRigidBodySimulation("simulation.mml");
    
    return 0;
}
```

---

## See Also

- [docs/MML_FILE_FORMATS.md](../docs/MML_FILE_FORMATS.md) - Complete format specifications
- [docs/core/Serializer.md](../docs/core/Serializer.md) - Serializer API documentation
- [docs/core/DataLoader.md](../docs/core/DataLoader.md) - DataLoader API documentation
- [tools/visualizers/](./visualizers/) - Visualization tools

---

*Last updated: February 2026*
*MML Version: 1.2+*
