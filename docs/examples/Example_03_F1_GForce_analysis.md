# Example: Formula 1 G-Force Analysis

**Where parametric curves meet real-world motorsport telemetry!**

<p align="center">
  <img src="../images/readme/examples/03_formula_1_sim/01_track_path.png" alt="F1 Track Layout" width="600">
</p>

## 📖 Overview

This example loads real Formula 1 telemetry data (position, speed, time) and uses MML's parametric curve and interpolation tools to compute the **G-forces** a driver experiences around the circuit — lateral forces in corners, longitudinal forces under braking and acceleration, and the combined total G-load.

> *"Under braking for Turn 1 at Monza, an F1 driver experiences over 5G — their body momentarily weighs five times what it does at rest."*

**Source:** `src/examples/03_formula_1_sim/`

**Build:** `cmake --build build --target Example03_F1GForces`

## 🏎️ Available Tracks

| Track | Driver | Length | Lap Time | Character |
|-------|--------|--------|----------|-----------|
| Silverstone (default) | Lewis Hamilton | ~5.8 km | ~89 s | High-speed sweeping corners |
| Monza (`--track=monza`) | Lando Norris | ~5.7 km | ~81 s | Long straights, heavy braking zones |

## 🎯 Analysis Pipeline

### Step 1: Load Telemetry Data

Reads a CSV file containing per-sample: time (s), distance (m), x/y position (m), and speed (km/h).

```
╔══════════════════════════════════════════════════════════╗
║         SILVERSTONE LAP TELEMETRY SUMMARY                ║
╠══════════════════════════════════════════════════════════╣
║  Samples:             871                                ║
║  Lap Time:          88.84 s                              ║
║  Track Length:      5790.6 m  (5.79 km)                  ║
║  Min Speed:           87.2 km/h                          ║
║  Max Speed:          325.4 km/h                          ║
║  Avg Speed:          234.7 km/h                          ║
╚══════════════════════════════════════════════════════════╝
```

### Step 2: Create Parametric Curve

Track coordinates (x, y) are fit with a `SplineInterpParametricCurve<2>`, parameterized by **distance along the track**. The spline smoothly interpolates between telemetry points, enabling curvature calculation at any location.

```cpp
// Create spline-interpolated parametric curve
// Parameter t goes from 0 to total distance (~5790 m)
MML::SplineInterpParametricCurve<2> trackCurve(minDist, maxDist, trackPoints, false);
```

The interpolation error is verified against original data points — typically sub-millimetre accuracy.

### Step 3: Speed Interpolation

Speed as a function of distance is created using `SplineInterpRealFunc`, enabling smooth evaluation of v(s) and its derivatives at any point along the track.

```cpp
MML::SplineInterpRealFunc speedInterp(distances, speeds_ms);
```

<p align="center">
  <img src="../images/readme/examples/03_formula_1_sim/03_speed_profile.png" alt="Speed Profile" width="600">
</p>

### Step 4: G-Force Calculation

At 500 equally-spaced sample points along the track, three G-force components are computed:

```
G-Force Summary:
┌─────────────────────────────────────────┐
│  Max Lateral G:        4.82 G          │ ← Cornering
│  Max Acceleration:     1.43 G          │ ← Throttle
│  Max Braking:          4.97 G          │ ← Brakes
│  Max Total G:          5.21 G          │ ← Combined
└─────────────────────────────────────────┘
```

Peak G-force locations identify the most demanding corners on the circuit.

### Step 5: Visualization

Three output files are generated:
- **Track layout** — 2D parametric curve visualization
- **G-force plot** — Lateral, longitudinal, and total G vs distance (3 curves)
- **Speed profile** — Speed (km/h) vs distance

<p align="center">
  <img src="../images/readme/examples/03_formula_1_sim/02_g_forces.png" alt="G-Forces vs Distance" width="600">
</p>

## ⚛️ Physics

### Lateral G-Force (Cornering)

When a car follows a curved path, it experiences centripetal acceleration proportional to the square of speed and the curvature of the path:

$$G_{\text{lat}} = \frac{v^2 \cdot \kappa}{g}$$

where:
- $v$ = speed in m/s
- $\kappa$ = curvature of the track (computed from the parametric curve)
- $g$ = 9.81 m/s²

### Longitudinal G-Force (Acceleration / Braking)

Changes in speed produce longitudinal forces. Using the chain rule to convert from distance-parameterized derivatives:

$$G_{\text{long}} = \frac{1}{g} \cdot \frac{dv}{dt} = \frac{1}{g} \cdot v \cdot \frac{dv}{ds}$$

Positive values indicate acceleration; negative values indicate braking.

### Total G-Force

The combined G-load magnitude:

$$G_{\text{total}} = \sqrt{G_{\text{lat}}^2 + G_{\text{long}}^2}$$

### Curvature from Parametric Curves

MML computes curvature using the standard formula for a 2D parametric curve $\mathbf{r}(s) = (x(s), y(s))$:

$$\kappa = \frac{|x'y'' - y'x''|}{(x'^2 + y'^2)^{3/2}}$$

This is computed by `Curves::getCurvature2D()` using automatic numerical differentiation of the spline.

## 🔧 MML Features Used

| Feature | Usage |
|---------|-------|
| `SplineInterpParametricCurve<2>` | Smooth 2D track representation from discrete points |
| `SplineInterpRealFunc` | Speed and G-force interpolation as functions of distance |
| `Curves::getCurvature2D()` | Curvature calculation along parametric curves |
| `Matrix<Real>` | Storage of track coordinate data |
| `Vector<Real>` | Distance, speed, and G-force arrays |
| `Visualizer` | Track layout, G-force, and speed profile plots |

## 📂 Data Format

Telemetry CSV files (`src/examples/03_formula_1_sim/data/`):

```csv
time_s,distance_m,x_m,y_m,speed_kmh
0.00,0.0,0.00,0.00,280.5
0.10,7.8,7.79,0.12,281.2
...
```

## 🚀 Running

```bash
# Default: Silverstone
./build/src/examples/Release/Example03_F1GForces

# Monza
./build/src/examples/Release/Example03_F1GForces --track=monza
```
