# Example 00: N-Body Gravity Simulation

**MML's flagship demonstration of numerical computing power!**

Simulate gravitational interactions between multiple bodies using Newton's law of universal gravitation. This example showcases MML's ODE solvers, vector mathematics, and visualization capabilities working together to solve a classic physics problem.

| Solar System Simulation | Star Cluster Collision |
|:-----------------------:|:----------------------:|
| ![Solar System](../images/_Example00_solar_system.png) | ![Star Clusters](../images/_Example00_star_clusters.png) |
| *Inner planets orbital mechanics* | *Two star clusters gravitational interaction* |

**Source:** [`src/examples/00_N_body_gravity/`](../../src/examples/00_N_body_gravity/)

---

## 📖 Table of Contents

- [Physics Background](#-physics-background)
- [Two Simulation Scenarios](#-two-simulation-scenarios)
- [Code Walkthrough](#-code-walkthrough)
- [MML Features Used](#-mml-features-used)
- [Running the Example](#-running-the-example)
- [Results & Visualization](#-results--visualization)
- [Conservation Laws](#-conservation-laws)

---

## 🔬 Physics Background

### Newton's Law of Universal Gravitation

Every particle in the universe attracts every other particle with a force proportional to the product of their masses and inversely proportional to the square of the distance between them:

$$\vec{F}_{ij} = -G \frac{m_i m_j}{|\vec{r}_{ij}|^2} \hat{r}_{ij}$$

Where:
- $G$ = gravitational constant (6.674 × 10⁻¹¹ m³ kg⁻¹ s⁻²)
- $m_i, m_j$ = masses of bodies i and j
- $\vec{r}_{ij}$ = position vector from body i to body j
- $\hat{r}_{ij}$ = unit vector in direction of $\vec{r}_{ij}$

### N-Body Problem

For N bodies, each body experiences gravitational attraction from all other N-1 bodies. The total force on body i is:

$$\vec{F}_i = \sum_{j \neq i} -G \frac{m_i m_j}{|\vec{r}_{ij}|^2} \hat{r}_{ij}$$

This creates a system of 6N coupled first-order ODEs (3 position + 3 velocity components per body).

### Conservation Laws

In an isolated gravitational system, three quantities are conserved:

| Quantity | Formula | Physical Meaning |
|----------|---------|------------------|
| **Total Energy** | $E = \sum_i \frac{1}{2}m_i v_i^2 - \sum_{i<j} \frac{G m_i m_j}{r_{ij}}$ | Kinetic + Potential |
| **Linear Momentum** | $\vec{P} = \sum_i m_i \vec{v}_i$ | Center of mass velocity |
| **Angular Momentum** | $\vec{L} = \sum_i m_i (\vec{r}_i \times \vec{v}_i)$ | Rotational motion |

These conservation laws provide excellent validation of numerical integration accuracy.

---

## 🌌 Two Simulation Scenarios

### Scenario 1: Solar System

**Real astronomical data!** Simulates the Sun with all 8 planets using actual masses and orbital distances.

| Body | Mass (Jupiter masses) | Orbital Distance (million km) |
|------|----------------------|------------------------------|
| Sun | 1047.35 | 0 (center) |
| Mercury | 1.66×10⁻⁷ | 57.91 |
| Venus | 2.45×10⁻⁶ | 108.21 |
| Earth | 3.00×10⁻⁶ | 149.60 |
| Mars | 3.23×10⁻⁷ | 227.92 |
| Jupiter | 1.0 | 778.57 |
| Saturn | 0.299 | 1433.53 |
| Uranus | 0.0457 | 2872.46 |
| Neptune | 0.0539 | 4495.06 |

**Units:** million km, Jupiter masses, years

### Scenario 2: Star Cluster Collision 🌟

**The most visually stunning scenario!** Two star clusters approach each other and undergo gravitational interaction, creating beautiful chaotic dynamics.

| Step 1: Initial Approach | Step 2: First Contact |
|:------------------------:|:---------------------:|
| ![Step 1](../images/_Example00_star_clusters_particle_step%201.png) | ![Step 2](../images/_Example00_star_clusters_particle_step%202.png) |
| *Two clusters approaching* | *Gravitational capture begins* |

| Step 3: Intermingling | Step 4: Post-Collision |
|:---------------------:|:----------------------:|
| ![Step 3](../images/_Example00_star_clusters_particle_step%203.png) | ![Step 4](../images/_Example00_star_clusters_particle_step%204.png) |
| *Stars exchange between clusters* | *New equilibrium forming* |

**Configuration:**
- Two separate clusters, each with 50 stars
- Central massive bodies in each cluster (mass = 5,000)
- Smaller stars (masses 1-10) in quasi-circular orbits
- Clusters on collision course with relative velocity
- Complexity: O(N²) = 10,000 force calculations per time step

This scenario demonstrates:
- **Gravitational scattering** — stars deflected by close encounters
- **Tidal disruption** — outer stars stripped from clusters
- **Energy exchange** — kinetic/potential energy redistribution
- **Chaotic dynamics** — sensitive dependence on initial conditions

---

## 💻 Code Walkthrough

### Self-Contained Physics Engine

The entire N-body simulation is contained in `NBodyGravity.h` (~670 lines) with no external dependencies beyond MML:

```cpp
// Key data structures
class GravityBodyState {
    Real _mass, _radius;
    Vec3Cart _position, _velocity;
    // ... accessors and constructors
};

class NBodyState {
    std::vector<GravityBodyState> _bodies;
    
    // Conservation law calculations
    Vec3Cart LinearMomentum() const;
    Vec3Cart AngularMomentumCM() const;
    Real TotalKineticEnergy() const;
    Real TotalPotentialEnergy() const;
};
```

### Solar System Configuration

Real astronomical data is built into the example:

```cpp
static NBodyGravitySimConfig Config2_Solar_system()
{
    // G in units: [(million km)³ / (Jupiter mass × year²)]
    double G = G_SI * M_jup * (year * year) / (million_km * million_km * million_km);
    
    NBodyGravitySimConfig config(G);
    
    // Add Sun at origin
    config.AddBody(M_sun / M_jup, Vec3Cart{0,0,0}, Vec3Cart{0,0,0}, "Yellow", 30);
    
    // Add planets with circular orbit velocities: v = sqrt(G * M_sun / r)
    for (int i = 0; i < SolarSystem::NumPlanets; i++) {
        const auto& p = SolarSystem::Planets[i];
        double dist = p.orbital_dist_km / 1e6;
        double v = std::sqrt(G * (M_sun / M_jup) / dist);
        
        config.AddBody(p.mass_jupiter, 
                      Vec3Cart{dist, 0, 0}, 
                      Vec3Cart{0, v, 10},
                      colors[i], radii[i]);
    }
    return config;
}
```

### Running the Simulation

```cpp
// Create configuration and simulator
NBodyGravitySimConfig config = NBodyGravityConfigGenerator::Config2_Solar_system();
NBodyGravitySimulator solver(config);

// Simulation parameters
Real duration = 2.0;        // 2 years
const int steps = 1001;
const Real dt = duration / steps;

// Integrate using Euler method
NBodyGravitySimulationResults result = solver.SolveEuler(dt, steps);

// Or use adaptive RK5 Cash-Karp for better accuracy
// NBodyGravitySimulationResults result = solver.SolveRK5(dt, steps);
```

### Visualization

Results can be visualized two ways:

```cpp
// Static trajectory plot (parametric curves)
result.VisualizeAsParamCurve("solar_system", Vector<int>{ 1, 2, 3, 4, 5, 6, 7, 8 });

// Animated particle simulation (3D with playback controls)
result.VisualizeAsParticleSimulation("solar_system_anim", 
                                      Vector<int>{ 0, 1, 2, 3, 4, 5, 6, 7, 8 }, dt);
```

---

## 🛠️ MML Features Used

This example showcases many MML components working together:

| Component | Usage |
|-----------|-------|
| `Vec3Cart` | 3D position and velocity vectors |
| `Vector<T>` | Dynamic arrays for ODE state |
| `IODESystem` | Interface for ODE system definition |
| `ODESystemSolver` | Euler and RK5 integrators |
| `Visualizer` | Real-time 3D visualization |
| `Serializer` | Save trajectories to files |

### ODE System Interface

The N-body problem is formulated as an ODE system:

```cpp
class NBodyGravityODESystem : public IODESystem {
    // State vector: [x1, y1, z1, x2, y2, z2, ..., vx1, vy1, vz1, vx2, ...]
    // Derivatives: [vx1, vy1, vz1, ..., ax1, ay1, az1, ...]
    
    void derivs(Real t, const Vector<Real>& state, Vector<Real>& derivs) {
        // Position derivatives = velocities
        // Velocity derivatives = accelerations from gravity
    }
};
```

---

## 🚀 Running the Example

### Build

```bash
cmake -B build
cmake --build build --target Example00_NBodyGravity
```

### Run

```bash
# Windows
.\build\src\examples\Release\Example00_NBodyGravity.exe

# Linux/macOS
./build/src/examples/Example00_NBodyGravity
```

### Output

The program runs all three scenarios and produces:
- Console output with energy conservation checks
- Trajectory files in `results/`:
  - `solar_system_trajectories.txt`
  - `five_body_trajectories.txt`
  - `many_body_trajectories.txt`
- Real-time visualizations (if visualizers are available)

---

## 📊 Results & Visualization

### Star Cluster Collision Sequence

The most spectacular visualization — watch two star clusters collide and interact:

| Overview | Particle Simulation |
|:--------:|:-------------------:|
| ![Cluster Overview](../images/_Example00_star_clusters_particle_overview.png) | ![Particle Sim](../images/_Example00_solar_system_particle_sim.png) |
| *Full cluster collision dynamics* | *Real-time particle visualization* |

### Platform-Specific Viewers

| Platform | Screenshot |
|----------|------------|
| **Windows** | ![Windows Particle Visualizer](../images/_Win%20-%20Particle%20visualizer%203D.png) |
| **Linux** | ![Linux Particle Visualizer](../images/_Linux%20-%20Particle%20visualizer%203D.png) |

### Sample Output

```
╔══════════════════════════════════════════════════════════════════════╗
║         MML N-BODY GRAVITATIONAL SIMULATOR                           ║
╚══════════════════════════════════════════════════════════════════════╝

======================================================================
  SCENARIO 1: Solar System Simulation
======================================================================

Solar System Configuration:
  Bodies: 9 (Sun + 8 planets)
  Gravitational constant G = 2.959122e-04
  Units: million km, Jupiter masses, years

Energy Conservation Check:
  Initial total energy: -1.234567e+08
  Final total energy:   -1.234520e+08
  Relative error: 0.0038%

✓ Solar system simulation complete!
```

---

## ⚖️ Conservation Laws

### Energy Conservation Analysis

The Euler method introduces numerical energy drift. For better conservation, use RK5:

| Method | Duration | Energy Drift |
|--------|----------|--------------|
| Euler | 2 years | ~0.004% |
| RK5 Cash-Karp | 2 years | ~10⁻¹⁰% |

### Momentum Conservation

```
Linear Momentum (should be constant):
  Initial: (0.55, -0.7, 0.9)
  Final:   (0.55, -0.7, 0.9)  ✓
```

---

## 📚 Further Reading

- **Numerical Recipes** - Chapter on N-body methods
- **Aarseth, S.J.** - *Gravitational N-Body Simulations* (Cambridge, 2003)
- **MML Documentation:**
  - [ODE Solvers](../algorithms/Differential_equations_solvers.md)
  - [Visualizers](../tools/Visualizers.md)

---

## 🔗 Related Examples

| Example | Description |
|---------|-------------|
| [Double Pendulum](Example_02_double_pendulum.md) | Another chaotic dynamics example |
| [Collision Simulator](Example_03_collision_simulator_2d.md) | Large-scale particle systems |
