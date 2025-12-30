# MML Examples — Physics Simulations That BLOW YOUR MIND! 🚀

**Five spectacular physics demonstrations** showcasing MML's capabilities.

## 🎯 Example Gallery

| # | Example | Physics | Bodies/Particles | Highlight |
|---|---------|---------|-----------------|-----------|
| 00 | [**N-Body Gravity**](Example_00_N_body_gravity.md) | Orbital mechanics | 200+ | Star cluster collision! |
| 01 | [**Projectile Launch**](Example_01_projectile_launch.md) | Ballistics | 1 | Air drag ruins everything |
| 02 | [**Double Pendulum**](Example_02_double_pendulum.md) | Chaos theory | 2 | Butterfly effect: 0.001°→47° |
| 03 | [**Collision Simulator**](Example_03_collision_simulator_2d.md) | Kinetic theory | 30,000+ | Shock wave propagation |
| 04 | [**Lorentz Transforms**](Example_04_Lorentz_transformations.md) | Special relativity | 2 twins | Time dilation is REAL |

## 🌌 Featured: Star Cluster Collision (Example 00)

Watch two 100-body star clusters undergo a hyperbolic flyby with gravitational deflection:

```
Configuration: Config3_StarClusterCollision
  - Two 101-body clusters (202 total)
  - Impact parameter: 250 units (not head-on!)
  - Approach speed: 5.0 units
  - Dramatic gravitational slingshot effect!
```

## 🔬 Featured: 30,000 Particle Shock Wave (Example 03)

100 high-energy particles explode into 29,900 cold particles:

```
Energy propagation (avg speed evolution):
  Step    Avg Speed
  ----    ---------
     0        5.16    ← Almost all cold
    80       23.45    ← Shock wave expanding
   160       35.67    ← Approaching equilibrium
```

## 📁 Source Locations

| Example | Source Directory |
|---------|-----------------|
| 00 | `src/examples/00_N_body_gravity/` |
| 01 | `src/examples/01_projectile_launch/` |
| 02 | `src/examples/02_double_pendulum/` |
| 03 | `src/examples/03_collision_simulator_2d/` |
| 04 | `src/examples/04_Lorentz_transformations/` |

## 🏃 Building All Examples

```bash
# Build all examples
cmake --build build --target examples

# Or build individually
cmake --build build --target Example00_NBodyGravity
cmake --build build --target Example01_ProjectileLaunch
cmake --build build --target Example02_DoublePendulum
cmake --build build --target Example03_CollisionSim2D
cmake --build build --target Example05_LorentzTransform
```

## ✨ Common Features

All examples demonstrate MML's core capabilities:
- **ODE Integration** — RK4, adaptive methods
- **Visualization** — Real-time particle and curve plotting
- **Serialization** — Export results for further analysis
- **Self-contained physics** — No external dependencies

## 📚 Start Here

New to MML? Start with **Example 01 (Projectile)** for a gentle introduction, then try **Example 02 (Double Pendulum)** for chaos, and finally tackle **Example 00 (N-Body)** for the full experience!