# Example 05: Rigid Body Collision Simulator

**Full 3D rigid body dynamics with quaternion rotations, elastic collisions, and conservation tracking.**

<p align="center">
  <img src="../images/readme/examples/05_rigid_body/01_rigid_body_start.png" alt="Rigid Body Simulation Start" width="500">
  <img src="../images/readme/examples/05_rigid_body/02_rigid_body.png" alt="Rigid Body Simulation" width="500">
</p>

## 📖 Overview

Three rigid bodies — two rectangular boxes and a sphere — tumble, spin, and collide inside a cubic container. The simulation tracks **13 degrees of freedom** per body (position, velocity, orientation quaternion, angular velocity) and verifies energy and momentum conservation at every step.

> *"The beauty of rigid body dynamics lies in the interplay between translational and rotational motion — a collision doesn't just change velocity, it transfers spin."*

**Source:** `src/examples/05_rigid_body/`

**Build:** `cmake --build build --target Example05_RigidBody`

## 🎯 Simulation Setup

### Bodies

| Body | Shape | Dimensions | Mass | Starting Position | Initial Velocity |
|------|-------|-----------|------|-------------------|------------------|
| Box 1 (red) | Parallelepiped | 2.0 × 1.0 × 0.6 m | 10 kg | (−2.0, 0.5, 0.0) | (3.0, 0.5, 0.2) m/s |
| Box 2 (blue) | Parallelepiped | 1.5 × 1.0 × 0.8 m | 8 kg | (2.0, −0.5, 0.3) | (−2.0, 0.3, −0.1) m/s |
| Sphere (green) | Sphere | r = 0.35 m | 5 kg | (−3.0, 3.5, −2.0) | (3.0, −6.0, 1.5) m/s |

All bodies also carry initial angular velocity and orientation (via `Quaternion::FromAxisAngle`).

### Environment

- **Container:** 10 m × 10 m × 10 m cubic box (walls at ±5 m)
- **Collisions:** Perfectly elastic (coefficient of restitution e = 1.0)
- **Time step:** dt = 0.001 s (1 ms)
- **Duration:** 40 seconds

### Configuration

```cpp
SimulationConfig config;
config.timeStep = 0.001;          // 1ms timestep
config.totalTime = 40.0;          // 40 second simulation
config.containerHalfSize = 5.0;   // 10m × 10m × 10m container
config.coeffRestitution = 1.0;    // Perfectly elastic
```

## 🧮 Conservation Tracking

The simulation continuously monitors three conserved quantities:

```
╔═══════════════════════════════════════════════════════════════════╗
║                        CONSERVATION CHECK                         ║
╠═══════════════════════════════════════════════════════════════════╣
║  KINETIC ENERGY:                                                  ║
║  Initial Total:   157.825000 J  (trans: 142.125, rot: 15.700)    ║
║  Final Total:     157.825000 J  (trans: 139.482, rot: 18.343)    ║
║  Energy Error:       0.000000 %  ✓ CONSERVED                     ║
╠═══════════════════════════════════════════════════════════════════╣
║  LINEAR MOMENTUM (vector sum changes on wall bounces):            ║
║  ANGULAR MOMENTUM (vector sum changes on wall bounces):           ║
╚═══════════════════════════════════════════════════════════════════╝
```

**Key insight:** Total kinetic energy (translational + rotational) is strictly conserved, but energy flows between types — collisions convert translational energy into rotational energy and vice versa.

Linear and angular momentum are conserved between body-body collisions, but wall bounces change the system total (the walls have effectively infinite mass).

## ⚛️ Physics

### Rigid Body State (13 DOF)

Each body carries a complete state vector:

| DOF | Quantity | Representation |
|-----|----------|----------------|
| 3 | Position | $\mathbf{r} = (x, y, z)$ |
| 3 | Velocity | $\mathbf{v} = (v_x, v_y, v_z)$ |
| 4 | Orientation | Quaternion $q = (w, x, y, z)$ with $\|q\| = 1$ |
| 3 | Angular velocity | $\boldsymbol{\omega} = (\omega_x, \omega_y, \omega_z)$ in body frame |

### Quaternion Orientation

Orientation is represented as a unit quaternion rather than Euler angles, avoiding gimbal lock:

$$q = w + xi + yj + zk, \quad \|q\| = 1$$

The rotation matrix $R$ converting from body frame to world frame is derived from $q$. Quaternion time evolution:

$$\dot{q} = \frac{1}{2} q \otimes \boldsymbol{\omega}$$

where $\boldsymbol{\omega}$ is treated as a pure quaternion $(0, \omega_x, \omega_y, \omega_z)$.

### Inertia Tensor

For a **box** with half-dimensions $(a, b, c)$ and mass $m$:

$$I = \begin{pmatrix} \frac{m}{3}(b^2 + c^2) & 0 & 0 \\ 0 & \frac{m}{3}(a^2 + c^2) & 0 \\ 0 & 0 & \frac{m}{3}(a^2 + b^2) \end{pmatrix}$$

For a **sphere** with radius $r$ and mass $m$:

$$I = \frac{2}{5} m r^2 \cdot \mathbf{I}_{3 \times 3}$$

### Kinetic Energy

$$E_{\text{total}} = \underbrace{\frac{1}{2} m |\mathbf{v}|^2}_{\text{translational}} + \underbrace{\frac{1}{2} \boldsymbol{\omega}^T I \boldsymbol{\omega}}_{\text{rotational}}$$

### Elastic Collision Response

At contact, impulse is applied to both bodies. For perfectly elastic collisions ($e = 1$), the impulse magnitude along the contact normal $\hat{n}$ is:

$$j = \frac{-(1 + e) \, v_{\text{rel}} \cdot \hat{n}}{\frac{1}{m_1} + \frac{1}{m_2} + \hat{n} \cdot (I_1^{-1}(\mathbf{r}_1 \times \hat{n})) \times \mathbf{r}_1 + \hat{n} \cdot (I_2^{-1}(\mathbf{r}_2 \times \hat{n})) \times \mathbf{r}_2}$$

This simultaneously updates both linear and angular velocities.

## 🔧 MML Features Used

| Feature | Usage |
|---------|-------|
| `Quaternion` | Orientation representation, `FromAxisAngle()` construction |
| `Vec3Cart` | Position, velocity, angular velocity vectors |
| `MatrixNM` | 3×3 inertia tensors in body frame |
| `Visualizer` | Launches WPF 3D rigid body visualizer |

### Self-Contained Headers

The rigid body simulation uses three local headers extracted from the library:

| Header | Purpose |
|--------|---------|
| `RigidBodyCore.h` | `RigidBody`, `RigidBodyBox`, `RigidBodySphere`, `RigidBodyState` |
| `RigidBodySimulator.h` | `RigidBodySimulator`, `SimulationConfig`, collision detection/response |
| `RigidBodySerializer.h` | `.mml` trajectory file output for the 3D visualizer |

## 📂 Output Files

The simulation produces trajectory files for 3D visualization:

| File | Contents |
|------|----------|
| `rigid_body_simulation.mml` | Combined trajectory (all bodies) |
| `rigid_body_box1.mml` | Box 1 trajectory only |
| `rigid_body_box2.mml` | Box 2 trajectory only |
| `rigid_body_sphere.mml` | Sphere trajectory only |

These `.mml` files can be loaded into the WPF-based 3D rigid body visualizer at `tools/visualizers/win/WPF/MML_RigidBodyMovement_Visualizer/`.

## 🚀 Running

```bash
# Build and run
cmake --build build --target Example05_RigidBody
./build/src/examples/Release/Example05_RigidBody
```

The simulation prints progress every second with current energy and drift percentage, then launches the 3D visualizer automatically.
