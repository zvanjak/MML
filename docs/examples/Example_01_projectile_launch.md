# Example 01: Projectile Launch with Air Resistance

**The classic physics problem, done RIGHT** - with real atmospheric models and drag!

<p align="center">
  <img src="../images/examples/projectile_vacuum_vs_air.png" alt="Vacuum vs Air Trajectories" width="600">
</p>

## 📖 Overview

This example demonstrates projectile motion simulation with MML, comparing:
- **Vacuum trajectories** (analytical and numerical)
- **Air resistance effects** (quadratic drag)
- **Optimal launch angles** (spoiler: NOT 45° with drag!)

**Source:** `src/examples/01_projectile_launch/`

## 🎯 Scenarios

### Scenario 1: Vacuum vs Air Resistance

**The fundamental question:** How much does air resistance matter?

```
Launch: v₀ = 100 m/s, θ = 45°

Vacuum (analytical):
  Range: 1019.37 m
  Max height: 254.84 m
  
With Air Resistance:
  Range: 432.15 m      ← 58% REDUCTION!
  Max height: 186.22 m
```

**Conclusion:** Air resistance is NOT negligible! At reasonable speeds, it can cut your range by more than half.

### Scenario 2: Optimal Launch Angle

In vacuum, the optimal angle is exactly 45°. But with air resistance...

| Angle | Range (m) |
|-------|-----------|
| 30°   | 176.8     |
| 35°   | 183.2     |
| **40°** | **185.7** ← Optimal! |
| 45°   | 183.5     |
| 50°   | 176.1     |

**The optimal angle drops to ~40°** because higher trajectories spend more time in the air, experiencing more drag.

### Scenario 3: Baseball Types

Different ball surfaces have different drag characteristics:
- Smooth ball: Lower drag coefficient
- Rough/stitched ball: Higher drag coefficient

## ⚛️ Physics

### Equations of Motion

**Without drag (vacuum):**
```
dx/dt = vₓ
dy/dt = vᵧ  
dvₓ/dt = 0
dvᵧ/dt = -g
```

**With quadratic drag:**
```
F_drag = -k|v|²v̂ = -k·v·v
    
dvₓ/dt = -k·|v|·vₓ / m
dvᵧ/dt = -g - k·|v|·vᵧ / m
```

Where:
- **k** = drag coefficient (depends on air density, cross-section, drag factor)
- **|v|** = speed = √(vₓ² + vᵧ²)
- **g** = 9.81 m/s² (gravitational acceleration)

### Analytical Solutions (Vacuum Only)

Range: $R = \frac{v_0^2 \sin(2\theta)}{g}$

Max height: $H = \frac{v_0^2 \sin^2(\theta)}{2g}$

Time of flight: $T = \frac{2v_0 \sin(\theta)}{g}$

## 🔧 MML Features Used

| Feature | Usage |
|---------|-------|
| `IODESystem` | Define projectile dynamics |
| `ODESystemSolver` | Integrate equations of motion |
| `Euler` / `RK4` | Numerical integration methods |
| `Visualizer` | Plot trajectories |
| `Serializer` | Export results |

## 📁 Key Files

```
src/examples/01_projectile_launch/
├── main.cpp              # Demo scenarios
├── Projectiles2D.h       # Self-contained projectile physics
└── CMakeLists.txt
```

## 🏃 Running

```bash
# Build
cmake --build build --target Example01_ProjectileLaunch

# Run
./build/src/examples/Release/Example01_ProjectileLaunch
```

## 📊 Sample Output

```
======================================================================
  SCENARIO 1: Vacuum vs Air Resistance
======================================================================

Launch Parameters:
  Initial velocity: 100 m/s
  Launch angle: 45 degrees
  Initial height: 0 m

Vacuum (analytical):
  Range: 1019.37 m
  Time of flight: 14.42 s
  Max height: 254.84 m

With Air Resistance (numerical):
  Range: 432.15 m
  Time of flight: 9.87 s
  Range reduction: 57.6%

✓ Vacuum vs Air comparison complete!
```

## 🎓 Learning Points

1. **Analytical vs Numerical:** Vacuum has closed-form solutions; drag requires numerical integration
2. **Physical intuition:** Air resistance is NOT negligible at everyday speeds
3. **ODE formulation:** Any physics problem can be cast as `dx/dt = f(x, t)`
4. **MML workflow:** Define system → Choose solver → Integrate → Visualize

## 📚 References

- Marion & Thornton, *Classical Dynamics* - Chapter on air resistance
- Halliday & Resnick, *Fundamentals of Physics* - Projectile motion
