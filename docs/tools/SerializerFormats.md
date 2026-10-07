# MML Serializer File Formats

This document describes the `.mml` presentation-export file formats produced by visualizer-oriented `Serializer` methods. All formats are plain text with consistent header structures.

These files are stable inputs for visualizers and examples, but they are not the
round-trip object-persistence formats. For `.mmlj` JSON and `.mmlb` binary
object persistence, see [SerializationPersistence.md](SerializationPersistence.md).

All type strings carry the `MML_` prefix and every header includes a `VERSION:` line immediately after the type identifier. The current version is **1**.

Type string constants and the version number are defined in `FormatType` namespace (see `SerializerBase.h`). Titles and legend entries are single-line display metadata; carriage returns and line feeds are written as spaces so accidental multiline input cannot corrupt the header structure.

---

## Real Function (1D)

**Type:** `MML_REAL_FUNCTION`  
**Functions:** `SaveRealFunc`, `SaveODESolutionComponentAsFunc`

```
MML_REAL_FUNCTION
VERSION: 1
sin(x)
x1: 0
x2: 6.28319
NumPoints: 100
0 0
0.0634665 0.0634239
0.126933 0.126592
...
```

---

## Real Function (Equally Spaced)

**Type:** `MML_REAL_FUNCTION_EQUALLY_SPACED`  
**Functions:** `SaveRealFuncEquallySpaced`

```
MML_REAL_FUNCTION_EQUALLY_SPACED
VERSION: 1
cos(x)
x1: 0
x2: 6.28319
NumPoints: 100
1
0.998027
0.992115
...
```

---

## Multi-Function (Multiple 1D Functions)

**Type:** `MML_MULTI_REAL_FUNCTION`  
**Functions:** `SaveRealMultiFunc`, `SaveODESolutionAsMultiFunc`

```
MML_MULTI_REAL_FUNCTION
VERSION: 1
Comparison
NumFuncs: 3
Legend: sin(x), cos(x), tan(x)
x1: 0
x2: 3.14159
NumPoints: 50
0 0 1 0
0.0641409 0.0640702 0.997945 0.0642925
...
```

---

## Parametric Curve 2D

**Type:** `MML_PARAMETRIC_CURVE_CARTESIAN_2D`  
**Functions:** `SaveParamCurve<2>`, `SaveAsParamCurve2D`, `SaveODESolAsParametricCurve2D`

```
MML_PARAMETRIC_CURVE_CARTESIAN_2D
VERSION: 1
Circle
t1: 0
t2: 6.28319
NumPoints: 100
0 1 0
0.0634665 0.997986 0.0634239
0.126933 0.991967 0.126592
...
```

Data columns: `t x(t) y(t)`

---

## Parametric Curve 3D

**Type:** `MML_PARAMETRIC_CURVE_CARTESIAN_3D`  
**Functions:** `SaveParamCurve<3>`, `SaveODESolAsParametricCurve3D`

```
MML_PARAMETRIC_CURVE_CARTESIAN_3D
VERSION: 1
Helix
t1: 0
t2: 12.5664
NumPoints: 200
0 1 0 0
0.0631653 0.998005 0.0631455 0.0100531
...
```

Data columns: `t x(t) y(t) z(t)`

---

## Scalar Function 2D (Surface)

**Type:** `MML_SCALAR_FUNCTION_CARTESIAN_2D`  
**Functions:** `SaveScalarFunc2DCartesian`

```
MML_SCALAR_FUNCTION_CARTESIAN_2D
VERSION: 1
z = x*y
x1: -5
x2: 5
NumPointsX: 50
y1: -5
y2: 5
NumPointsY: 50
-5 -5 25
-5 -4.79592 23.9796
...
```

Data columns: `x y f(x,y)`

---

## Scalar Function 3D (Volumetric)

**Type:** `MML_SCALAR_FUNCTION_CARTESIAN_3D`  
**Functions:** `SaveScalarFunc3DCartesian`

```
MML_SCALAR_FUNCTION_CARTESIAN_3D
VERSION: 1
w = x*y*z
x1: -2
x2: 2
NumPointsX: 10
y1: -2
y2: 2
NumPointsY: 10
z1: -2
z2: 2
NumPointsZ: 10
-2 -2 -2 -8
-2 -2 -1.55556 -6.22222
...
```

Data columns: `x y z f(x,y,z)`

---

## Vector Field 2D

**Type:** `MML_VECTOR_FIELD_2D_CARTESIAN`  
**Functions:** `SaveVectorFunc2D`, `SaveVectorFunc2DCartesian`

```
MML_VECTOR_FIELD_2D_CARTESIAN
VERSION: 1
Rotation field
-5 -5 5 -5
-5 -4 4 -5
-5 -3 3 -5
...
```

Data columns: `x y Fx(x,y) Fy(x,y)`

---

## Vector Field 3D

**Type:** `MML_VECTOR_FIELD_3D_CARTESIAN`  
**Functions:** `SaveVectorFunc3D`, `SaveVectorFunc3DCartesian`

```
MML_VECTOR_FIELD_3D_CARTESIAN
VERSION: 1
Gravity field
-2 -2 -2 0.096225 0.096225 0.096225
-2 -2 -1 0.111111 0.111111 0.0555556
...
```

Data columns: `x y z Fx Fy Fz`

---

## Particle Simulation 2D

**Type:** `MML_PARTICLE_SIMULATION_DATA_2D`  
**Functions:** `SaveParticleSimulation2D`

```
MML_PARTICLE_SIMULATION_DATA_2D
VERSION: 1
Width: 800
Height: 600
NumBalls: 3
Ball_1 red 5
Ball_2 blue 5
Ball_3 green 5
NumSteps: 100
Step 0 0
0 100.5 200.3
1 150.2 300.1
2 250.0 150.5
Step 1 0.016
0 101.2 201.5
1 151.0 301.2
2 251.3 151.8
...
```

Each step contains one row per ball: `ballIndex x y`

---

## Particle Simulation 3D

**Type:** `MML_PARTICLE_SIMULATION_DATA_3D`  
**Functions:** `SaveParticleSimulation3D`

```
MML_PARTICLE_SIMULATION_DATA_3D
VERSION: 1
Width: 400
Height: 400
Depth: 400
NumBalls: 3
Ball_1 red 10
Ball_2 blue 10
Ball_3 green 10
NumSteps: 100
Step 0 0
0 100 200 150
1 150 300 200
2 250 150 100
Step 1 0.016
0 101 201 151
1 151 301 201
2 251 151 101
...
```

Each step contains one row per ball: `ballIndex x y z`

---

## Vector Field (Spherical)

**Type:** `MML_VECTOR_FIELD_SPHERICAL`  
**Functions:** `SaveVectorFuncSpherical`

```
MML_VECTOR_FIELD_SPHERICAL
VERSION: 1
<title>
<r> <theta> <phi> <Fr> <Ftheta> <Fphi>
...
```

Data columns: `r θ φ F_r F_θ F_φ`

---

## Field Lines 2D

**Type:** `MML_FIELD_LINES_2D`  
**Functions:** `VisualizeFieldLines2D` (via Visualizer)

```
MML_FIELD_LINES_2D
VERSION: 1
<title>
<line data>
...
```

---

## Field Lines 3D

**Type:** `MML_FIELD_LINES_3D`  
**Functions:** `VisualizeFieldLines3D` (via Visualizer)

```
MML_FIELD_LINES_3D
VERSION: 1
<title>
<line data>
...
```

---

## Complete Type String Reference

All format type identifiers are defined as `constexpr const char*` constants in the `FormatType` namespace (`SerializerBase.h`):

| Constant | String Value | Version |
|----------|-------------|---------|
| `FormatType::REAL_FUNCTION` | `MML_REAL_FUNCTION` | 1 |
| `FormatType::REAL_FUNCTION_EQUALLY_SPACED` | `MML_REAL_FUNCTION_EQUALLY_SPACED` | 1 |
| `FormatType::MULTI_REAL_FUNCTION` | `MML_MULTI_REAL_FUNCTION` | 1 |
| `FormatType::PARAMETRIC_CURVE_CARTESIAN_2D` | `MML_PARAMETRIC_CURVE_CARTESIAN_2D` | 1 |
| `FormatType::PARAMETRIC_CURVE_CARTESIAN_3D` | `MML_PARAMETRIC_CURVE_CARTESIAN_3D` | 1 |
| `FormatType::PARAMETRIC_SURFACE_CARTESIAN` | `MML_PARAMETRIC_SURFACE_CARTESIAN` | 1 |
| `FormatType::SCALAR_FUNCTION_CARTESIAN_2D` | `MML_SCALAR_FUNCTION_CARTESIAN_2D` | 1 |
| `FormatType::SCALAR_FUNCTION_CARTESIAN_3D` | `MML_SCALAR_FUNCTION_CARTESIAN_3D` | 1 |
| `FormatType::VECTOR_FIELD_2D_CARTESIAN` | `MML_VECTOR_FIELD_2D_CARTESIAN` | 1 |
| `FormatType::VECTOR_FIELD_3D_CARTESIAN` | `MML_VECTOR_FIELD_3D_CARTESIAN` | 1 |
| `FormatType::VECTOR_FIELD_SPHERICAL` | `MML_VECTOR_FIELD_SPHERICAL` | 1 |
| `FormatType::FIELD_LINES_2D` | `MML_FIELD_LINES_2D` | 1 |
| `FormatType::FIELD_LINES_3D` | `MML_FIELD_LINES_3D` | 1 |
| `FormatType::PARTICLE_SIMULATION_DATA_2D` | `MML_PARTICLE_SIMULATION_DATA_2D` | 1 |
| `FormatType::PARTICLE_SIMULATION_DATA_3D` | `MML_PARTICLE_SIMULATION_DATA_3D` | 1 |

The version number is available as `FormatType::CURRENT_VERSION` (currently `1`).