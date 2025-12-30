# MML Serializer File Formats

This document describes the file formats produced by `Serializer` class methods. All formats are plain text with consistent header structures.

---

## Real Function (1D)

**Type:** `REAL_FUNCTION`  
**Functions:** `SaveRealFunc`, `SaveODESolutionComponentAsFunc`

```
REAL_FUNCTION
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

**Type:** `REAL_FUNCTION_EQUALLY_SPACED`  
**Functions:** `SaveRealFuncEquallySpaced`

```
REAL_FUNCTION_EQUALLY_SPACED
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

**Type:** `REAL_FUNCTION_MULTI`  
**Functions:** `SaveRealMultiFunc`, `SaveODESolutionAsMultiFunc`

```
REAL_FUNCTION_MULTI
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

**Type:** `PARAMETRIC_CURVE_CARTESIAN_2D`  
**Functions:** `SaveParamCurve<2>`, `SaveAsParamCurve2D`, `SaveODESolAsParametricCurve2D`

```
PARAMETRIC_CURVE_CARTESIAN_2D
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

**Type:** `PARAMETRIC_CURVE_CARTESIAN_3D`  
**Functions:** `SaveParamCurve<3>`, `SaveODESolAsParametricCurve3D`

```
PARAMETRIC_CURVE_CARTESIAN_3D
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

**Type:** `SCALAR_FUNCTION_CARTESIAN_2D`  
**Functions:** `SaveScalarFunc2DCartesian`

```
SCALAR_FUNCTION_CARTESIAN_2D
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

**Type:** `SCALAR_FUNCTION_CARTESIAN_3D`  
**Functions:** `SaveScalarFunc3DCartesian`

```
SCALAR_FUNCTION_CARTESIAN_3D
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

**Type:** `VECTOR_FIELD_2D_CARTESIAN`  
**Functions:** `SaveVectorFunc2D`, `SaveVectorFunc2DCartesian`

```
VECTOR_FIELD_2D_CARTESIAN
Rotation field
-5 -5 5 -5
-5 -4 4 -5
-5 -3 3 -5
...
```

Data columns: `x y Fx(x,y) Fy(x,y)`

---

## Vector Field 3D

**Type:** `VECTOR_FIELD_3D_CARTESIAN`  
**Functions:** `SaveVectorFunc3D`, `SaveVectorFunc3DCartesian`

```
VECTOR_FIELD_3D_CARTESIAN
Gravity field
-2 -2 -2 0.096225 0.096225 0.096225
-2 -2 -1 0.111111 0.111111 0.0555556
...
```

Data columns: `x y z Fx Fy Fz`

---

## Particle Simulation 2D

**Type:** `PARTICLE_SIMULATION_DATA_2D`  
**Functions:** `SaveParticleSimulation2D`

```
PARTICLE_SIMULATION_DATA_2D
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

**Type:** `PARTICLE_SIMULATION_DATA_3D`  
**Functions:** `SaveParticleSimulation3D`

```
PARTICLE_SIMULATION_DATA_3D
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