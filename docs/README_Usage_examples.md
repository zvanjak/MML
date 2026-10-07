# 🧪 MML Usage Examples — Physics Simulations You Can Run

> **Self-contained, production-ready examples** demonstrating MML's capabilities in real-world physics simulations. Each example includes all physics code in-directory — no external dependencies.

> Back to the [main README](../README.md) • See also [Code examples](README_Code_examples.md) and the [Visualization suite](README_Visualization_suite.md).

---

## 🌌 [Example 00: N-Body Gravity](examples/Example_00_N_body_gravity.md) — *Flagship Example*

**Solar System & Star Cluster Simulations** — Newton's law of universal gravitation with **7 integrators** (Euler, RK4, Verlet, Leapfrog, RK5, DP5, DP8). Symplectic integrators for long-term orbital stability, adaptive methods for highest accuracy. Self-contained physics engine (~870 lines).

<table>
<tr>
<td align="center" width="33%">

![Solar System](images/readme/examples/00_N_body_gravity/_Example00_solar_system.png)

*Solar system orbital mechanics*

</td>
<td align="center" width="33%">

![Particle Sim](images/readme/examples/00_N_body_gravity/_Example00_solar_system_particle_sim.png)

*Real-time particle visualization*

</td>
<td align="center" width="33%">

![Cluster Overview](images/readme/examples/00_N_body_gravity/_Example00_star_clusters_particle_overview.png)

*Star cluster collision overview*

</td>
</tr>
<tr>
<td align="center">

![Cluster Step 1](images/readme/examples/00_N_body_gravity/_Example00_star_clusters_particle_step%201.png)

*Cluster approach*

</td>
<td align="center">

![Cluster Step 3](images/readme/examples/00_N_body_gravity/_Example00_star_clusters_particle_step%203.png)

*Gravitational interaction*

</td>
<td align="center">

![Cluster Trajectories](images/readme/examples/00_N_body_gravity/_Example00_star_clusters_trajectories_visualization.png)

*Full trajectory visualization*

</td>
</tr>
</table>

```cpp
NBodyGravitySimConfig config = NBodyGravityConfigGenerator::Config1_Solar_system();
NBodyGravitySimulator solver(config);

auto results_verlet = solver.SolveVerlet(0.01, 365000);       // 10 years, symplectic
auto results_dp8 = solver.SolveDP8(10.0, 1e-12, 0.01, 0.001); // 10 years, adaptive DP8
```

---

## 🎾 [Example 02: Double Pendulum](examples/Example_02_double_pendulum.md) — *Chaos Theory*

Deterministic chaos in action — infinitely small differences in initial conditions lead to completely different outcomes.

<table>
<tr>
<td align="center" width="33%">

![Trajectory](images/readme/examples/02_double_pendulum/01_trajectory_for%20both%20angles.png)

*Angle trajectories*

</td>
<td align="center" width="33%">

![Phase Space](images/readme/examples/02_double_pendulum/02_phase_space.png)

*Phase space portrait*

</td>
<td align="center" width="33%">

![Butterfly Effect](images/readme/examples/02_double_pendulum/03_butterfly_effect.png)

*Butterfly effect divergence*

</td>
</tr>
</table>

---

## 🏎️ [Example 03: Formula 1 G-Force Analysis](examples/Example_03_F1_GForce_analysis.md) — *Parametric Curves*

Real F1 telemetry data (Silverstone, Monza) analyzed using MML's parametric curve and curvature calculations. Lateral G = v²κ/g, longitudinal G = (1/g)·dv/dt.

<table>
<tr>
<td align="center" width="33%">

![Track Path](images/readme/examples/03_formula_1_sim/01_track_path.png)

*Track layout from telemetry*

</td>
<td align="center" width="33%">

![G-Forces](images/readme/examples/03_formula_1_sim/02_g_forces.png)

*G-force profile around lap*

</td>
<td align="center" width="33%">

![Speed Profile](images/readme/examples/03_formula_1_sim/03_speed_profile.png)

*Speed profile analysis*

</td>
</tr>
</table>

---

## 💥 [Example 04: 2D Collision Simulator](examples/Example_04_collision_simulator_2d.md) — *Kinetic Theory*

**30,000+ particles** with exact elastic collision physics, spatial partitioning for O(N) performance, and multi-threaded execution. Watch shock waves propagate!

<table>
<tr>
<td align="center" width="33%">

![Shock 1](images/readme/examples/04_collision_simulator_2d/01_shock_wave.png)

*Initial shock front*

</td>
<td align="center" width="33%">

![Shock 2](images/readme/examples/04_collision_simulator_2d/02_shock_wave.png)

*Wave propagation*

</td>
<td align="center" width="33%">

![Shock 3](images/readme/examples/04_collision_simulator_2d/03_shock_wave.png)

*Shock wave dispersion*

</td>
</tr>
</table>

---

## 📦 [Example 05: Rigid Body Collisions](examples/Example_05_rigid_body.md) — *3D Dynamics*

Two parallelepipeds and a sphere in a cubic container with elastic collisions, full rotational dynamics using quaternions and inertia tensors.

<table>
<tr>
<td align="center" width="50%">

![Start](images/readme/examples/05_rigid_body/01_rigid_body_start.png)

*Initial configuration*

</td>
<td align="center" width="50%">

![Collision](images/readme/examples/05_rigid_body/02_rigid_body.png)

*Mid-collision dynamics*

</td>
</tr>
</table>

---

## More Examples

| # | Example | Description |
|---|---------|-------------|
| 01 | [**Projectile Launch**](examples/Example_01_projectile_launch.md) | Ballistic trajectory with air resistance — drag models, vacuum vs air comparison |
| 06 | [**Lorentz Transformations**](examples/Example_06_Lorentz_transformations.md) | Special relativity — acceleration, worldlines, proper time, and the Twin Paradox |

---

## 🚀 Try It Now

```bash
cmake -B build && cmake --build build

# Run the flagship N-Body simulation
./build/src/examples/Release/Example00_NBodyGravity      # Windows
./build/src/examples/Example00_NBodyGravity              # Linux
```

All examples produce visualization output viewable with the included Qt-based viewers. See the [Visualization suite](README_Visualization_suite.md) for the full gallery.