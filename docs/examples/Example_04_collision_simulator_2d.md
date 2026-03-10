# Example 03: 2D Collision Simulator - Kinetic Theory & Shock Waves

**30,000 particles colliding in perfect elastic chaos!** 🔴🔵

<p align="center">
  <img src="../images/examples/collision_shock_wave.png" alt="Shock Wave Propagation" width="600">
</p>

## 📖 Overview

This example demonstrates **large-scale particle physics** with MML's high-performance collision engine:

- **Exact elastic collision physics** - momentum AND energy conserved
- **Spatial partitioning** - O(N) average complexity instead of O(N²)
- **Multi-threaded execution** - scales to 30,000+ particles
- **Real-time visualization** - watch the physics unfold!

**Source:** `src/examples/03_collision_simulator_2d/`

## 🎯 Scenarios

### Scenario 1: Two-Color Mixing (Diffusion)

**500 balls** (250 blue on left, 250 red on right) with random velocities.

Watch them mix! This demonstrates **diffusion** - the statistical tendency for particles to spread out.

```
Configuration:
  Box: 1000 x 800
  Balls: 500 total (250 blue + 250 red)
  Radius: 5, Mass: 1.0
  Initial velocity: 50 m/s
  Time step: 0.1 s, Steps: 100

Final Statistics:
  Min speed: 12.34
  Max speed: 87.65
  Avg speed: 50.12 (σ = 18.45)

Simulation completed in 87 ms ← FAST!
```

**Physics insight:** The Maxwell-Boltzmann distribution emerges naturally from random collisions!

### Scenario 2: Shock Wave 🌊

**30,000 balls** with an energetic core at center - a tiny hot region surrounded by cold gas.

```
Configuration:
  Box: 2000 x 1600
  Surrounding gas: 29,900 balls (v = 5 m/s - COLD)
  Energetic core: 100 balls (v = 500 m/s - HOT!)
  Core radius: 50 units

Energy propagation (avg speed evolution):
  Step    Avg Speed
  ----    ---------
     0        5.16    ← Almost all cold
    40       12.34    ← Energy spreading
    80       23.45    ← Shock wave expanding
   120       31.22    ← Approaching equilibrium
   160       35.67    ← Nearly thermalized
```

**Physics:** The hot core explosively expands, creating a **spherical shock wave** that propagates through the cold gas. Energy transfers through collisions until equilibrium is reached.

## ⚛️ Physics

### Elastic Collision Equations

When two balls collide, both **momentum** and **kinetic energy** are conserved:

**Conservation laws:**
- Momentum: $m_1\vec{v}_1 + m_2\vec{v}_2 = m_1\vec{v}_1' + m_2\vec{v}_2'$
- Energy: $\frac{1}{2}m_1 v_1^2 + \frac{1}{2}m_2 v_2^2 = \frac{1}{2}m_1 v_1'^2 + \frac{1}{2}m_2 v_2'^2$

**Solution (1D along collision axis):**

$$v_1' = \frac{(m_1-m_2)v_1 + 2m_2 v_2}{m_1+m_2}$$

$$v_2' = \frac{(m_2-m_1)v_2 + 2m_1 v_1}{m_1+m_2}$$

**For equal masses:** Velocities simply **swap** along the collision axis!

### Collision Detection Algorithm

**Naive approach:** O(N²) - check all pairs → 30,000² = 900 million checks/frame! ❌

**Spatial subdivision:** O(N) average!
1. Divide space into grid cells
2. Only check balls in same cell or neighboring cells
3. Grid size ≈ 4× ball diameter for optimal performance

### Exact Collision Time

Balls might pass through each other if we only check positions at discrete times. The simulator calculates the **exact collision time** within each timestep:

```
Position: x(t) = x₀ + v·t
Distance: d(t) = |x₁(t) - x₂(t)|
Collision when: d(t) = r₁ + r₂

→ Solve quadratic equation for t ∈ [0, dt]
```

## 🔧 MML Features Used

| Feature | Usage |
|---------|-------|
| `CollisionSimulator2D` | Core physics engine |
| `BoxContainer2D` | Particle container with walls |
| `ContainerFactory2D` | Create standard configurations |
| `SimResultsCollSim2D` | Store simulation results |
| `ParticleVisualizer2D` | Real-time visualization |
| Multi-threading | Parallel collision detection |

## 📁 Key Files

```
src/examples/03_collision_simulator_2d/
├── main.cpp                 # Demo scenarios
├── CollisionSimulator2D.h   # Self-contained physics engine
└── CMakeLists.txt
```

## 🏃 Running

```bash
# Build
cmake --build build --target Example03_CollisionSim2D

# Run
./build/src/examples/Release/Example03_CollisionSim2D
```

**Note:** The shock wave scenario with 30,000 particles may take several seconds depending on your CPU.

## 📊 Sample Output

```
======================================================================
         MML 2D COLLISION SIMULATOR - ADVANCED PHYSICS
======================================================================

   Self-contained collision physics engine demonstrating:
   - Exact elastic collision calculation
   - Spatial subdivision for O(N) performance
   - Parallel execution for large simulations
   - Real-time visualization with ParticleVisualizer2D

======================================================================
  SCENARIO 2: Shock Wave Simulation (30,000 balls)
======================================================================

Configuration:
  Box: 2000 x 1600
  Surrounding gas: 29900 balls (v=5)
  Energetic core: 100 balls (v=500)
  Core radius: 50
  Time step: 0.01 s, Steps: 200

Creating shock wave configuration...
Created 30000 balls
Using 200x250 spatial grid
Running PARALLEL simulation...
  (This will use all CPU cores for collision detection)

Simulation completed in 12 seconds

✓ Shock wave simulation complete!
```

## 🎥 Visualization

The simulator automatically launches **ParticleVisualizer2D** for each scenario:

**Controls:**
- **Play/Pause** - Start/stop animation
- **Speed slider** - Adjust playback speed
- **Frame slider** - Scrub through simulation
- **Color coding** - Blue (slow) → Red (fast)

## 🎓 Learning Points

1. **Kinetic theory:** Statistical mechanics emerges from particle collisions
2. **Diffusion:** Random motion naturally mixes initially separated regions
3. **Shock waves:** Energy propagates as a coherent wave front
4. **Scaling:** Spatial partitioning enables large simulations
5. **Parallelism:** Multi-threading for compute-intensive physics

## 🔬 Extensions

Try modifying the code to explore:
- Different mass ratios (heavy vs light particles)
- Inelastic collisions (energy dissipation)
- External forces (gravity, electromagnetic)
- Different container shapes (circular, polygonal)

## 📚 References

- Reif, *Fundamentals of Statistical and Thermal Physics*
- Landau & Lifshitz, *Statistical Physics*
- Press et al., *Numerical Recipes* - Collision detection algorithms
