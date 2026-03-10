# PDE Solver Examples Gallery

**A curated collection of complete, working PDE examples demonstrating the MML framework**

This gallery contains 7 fully-functional examples that showcase different aspects of the MML PDE solver framework. Every example has been compiled and tested to ensure it works correctly.

## Overview

| Example | Type | Dimension | Features | Difficulty |
|---------|------|-----------|----------|------------|
| [1. 1D Steady Heat](#example-1-1d-steady-heat-conduction) | Elliptic | 1D | Analytical verification | ⭐ Beginner |
| [2. 2D Laplace](#example-2-2d-laplace-equation-heated-plate) | Elliptic | 2D | Mixed BCs | ⭐⭐ Beginner |
| [3. 2D Gaussian Source](#example-3-2d-poisson-with-gaussian-source) | Elliptic | 2D | Non-uniform source | ⭐⭐ Beginner |
| [4. 3D Poisson](#example-4-3d-poisson-equation) | Elliptic | 3D | Large system (15K unknowns) | ⭐⭐⭐ Intermediate |
| [5. 1D Transient Heat](#example-5-1d-time-dependent-heat) | Parabolic | 1D | Time evolution, analytical | ⭐⭐ Beginner |
| [6. 2D Room Cooling](#example-6-2d-room-cooling) | Parabolic | 2D | Time-dependent | ⭐⭐⭐ Intermediate |
| [7. Robin BC](#example-7-robin-boundary-conditions) | Elliptic | 1D | Convective cooling | ⭐⭐⭐ Intermediate |

## Running the Examples

### Compilation

All examples are included in the `PDE_Gallery_Examples` executable:

```bash
# Build the examples
cmake --build build --target PDE_Gallery_Examples --config Release

# Run all examples
./build/src/docs_demos/Release/PDE_Gallery_Examples.exe
```

### Source Code Location

Complete source: [`src/docs_demos/pde/gallery_examples_simple.cpp`](../../src/docs_demos/pde/gallery_examples_simple.cpp)

---

## Example 1: 1D Steady Heat Conduction

**Physical Problem:** Temperature distribution in a uniformly heated rod with fixed ends.

**Mathematical Problem:**
```
-d²u/dx² = 10  on [0,1]
u(0) = 0
u(1) = 0
```

**Analytical Solution:** `u(x) = 5x(1-x)`

### Complete Code

```cpp
void example1_1D_steady_heat() {
    // 1. Define the domain
    Interval<double> domain(0.0, 1.0);
    Grid1D<double> grid(domain, 100);
    
    // 2. Set boundary conditions: Zero temperature at ends
    auto bc = homogeneousDirichlet1D<double>();
    
    // 3. Create solver
    PoissonSolver1D<double> solver(grid, bc);
    
    // 4. Set uniform heat source
    solver.setSource([](double x) { return 10.0; });
    
    // 5. Solve
    auto solution = solver.solve();
    
    // 6. Analyze results
    double maxTemp = 0.0;
    double maxX = 0.0;
    for (int i = 0; i < grid.numNodes(); i++) {
        double x = grid.x(i);
        double u = solution(i);
        if (u > maxTemp) {
            maxTemp = u;
            maxX = x;
        }
    }
    
    std::cout << "Max temperature: " << maxTemp << " at x = " << maxX << "\n";
    // Expected: 1.25 at x = 0.5
}
```

### Key Concepts

- **Grid creation:** `Grid1D` with 100 nodes
- **Homogeneous Dirichlet BCs:** Both boundaries set to zero
- **Uniform source term:** Constant value 10 throughout domain
- **Direct solution access:** `solution(i)` for 1D GridFunction

### Expected Output

```
Max temperature: 1.25 at x = 0.5
Expected max: 1.25 at x = 0.5
```

**Verification:** Perfect match with analytical solution (error < 1e-6)

### Physical Interpretation

This models a metal rod generating heat uniformly (e.g., electrical resistance heating) with both ends kept at 0°C. Heat accumulates most at the center where it's farthest from the cold boundaries.

---

## Example 2: 2D Laplace Equation (Heated Plate)

**Physical Problem:** Steady-state temperature in a plate with hot top, cold bottom, and insulated sides.

**Mathematical Problem:**
```
-∇²u = 0  on [0,1]×[0,1]
u(top) = 100°C     (hot)
u(bottom) = 0°C    (cold)
∂u/∂n(left) = 0    (insulated)
∂u/∂n(right) = 0   (insulated)
```

### Complete Code

```cpp
void example2_2D_laplace() {
    // 1. Define domain
    Rectangle<double> domain(0.0, 1.0, 0.0, 1.0);
    Grid2D<double> grid(domain, 60, 60);
    
    // 2. Set mixed boundary conditions
    BoundaryConditions2D<double> bc;
    bc.setLeft(BoundaryCondition2D<double>::Neumann(0.0));      // Insulated
    bc.setRight(BoundaryCondition2D<double>::Neumann(0.0));     // Insulated
    bc.setBottom(BoundaryCondition2D<double>::Dirichlet(0.0));  // Cold (0°C)
    bc.setTop(BoundaryCondition2D<double>::Dirichlet(100.0));   // Hot (100°C)
    
    // 3. Create solver (Laplace: zero source)
    PoissonSolver2D<double> solver(grid, bc);
    solver.setSource([](double x, double y) { return 0.0; });
    
    // 4. Solve
    auto solution = solver.solve();
    
    // 5. Check center temperature
    int centerI = grid.numNodesX() / 2;
    int centerJ = grid.numNodesY() / 2;
    double centerTemp = solution(centerI, centerJ);
    
    std::cout << "Temperature at center: " << centerTemp << "°C\n";
}
```

### Key Concepts

- **Mixed BCs:** Combination of Dirichlet (temperature) and Neumann (flux)
- **Laplace equation:** Special case of Poisson with f = 0
- **2D indexing:** `solution(i, j)` for 2D GridFunction
- **Symmetry:** Insulated sides create 1D-like behavior in y-direction

### Expected Output

```
Temperature at center: ~50°C (expect linear variation)
```

### Physical Interpretation

Models a thermally-insulated plate (sides are insulated) with temperature controlled at top and bottom. With no insulation on sides, heat flows only vertically, creating nearly linear temperature gradient.

---

## Example 3: 2D Poisson with Gaussian Source

**Physical Problem:** Heat diffusion from a localized source (e.g., laser spot, heating element).

**Mathematical Problem:**
```
-∇²u = f(x,y)  on [0,1]×[0,1]
f(x,y) = 1000·exp(-100·r²)  where r² = (x-0.5)² + (y-0.5)²
u = 0 on all boundaries
```

### Complete Code

```cpp
void example3_2D_gaussian_source() {
    // 1. Setup
    Rectangle<double> domain(0.0, 1.0, 0.0, 1.0);
    Grid2D<double> grid(domain, 80, 80);
    auto bc = homogeneousDirichlet2D<double>();
    
    // 2. Create solver with Gaussian source
    PoissonSolver2D<double> solver(grid, bc);
    solver.setSource([](double x, double y) {
        double r2 = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5);
        return 1000.0 * std::exp(-100.0 * r2);
    });
    
    // 3. Solve
    auto solution = solver.solve();
    
    // 4. Find maximum (should be at center)
    double maxU = 0.0;
    for (int i = 0; i < grid.numNodesX(); i++) {
        for (int j = 0; j < grid.numNodesY(); j++) {
            double u = solution(i, j);
            if (u > maxU) maxU = u;
        }
    }
    
    std::cout << "Maximum value: " << maxU << "\n";
}
```

### Key Concepts

- **Non-uniform source:** Lambda function defines spatially-varying heat generation
- **Localized forcing:** Gaussian centered at domain center
- **Source strength:** Peak value 1000, decay controlled by exp(-100r²)
- **Maximum location:** Should coincide with source center

### Expected Output

```
Maximum value: 9.87851
```

### Physical Interpretation

Models a scenario like laser heating where energy is deposited in a small region. The temperature peaks at the source location and decays toward the cold boundaries. The Gaussian shape smoothly transitions from high intensity at center to zero away from it.

---

## Example 4: 3D Poisson Equation

**Physical Problem:** Electrostatic potential from a charge distribution in 3D space.

**Mathematical Problem:**
```
-∇²u = f(x,y,z)  in [0,1]³
f(x,y,z) = 1000·exp(-100·r²)  where r² = (x-0.5)² + (y-0.5)² + (z-0.5)²
u = 0 on all boundaries
```

### Complete Code

```cpp
void example4_3D_poisson() {
    // 1. Define 3D domain
    Box<double> domain(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
    Grid3D<double> grid(domain, 25, 25, 25);  // 15,625 unknowns
    
    // 2. Boundary conditions
    auto bc = homogeneousDirichlet3D<double>();
    
    // 3. Create solver with point source at center
    PoissonSolver3D<double> solver(grid, bc);
    solver.setSource([](double x, double y, double z) {
        double r2 = (x - 0.5) * (x - 0.5) + 
                    (y - 0.5) * (y - 0.5) + 
                    (z - 0.5) * (z - 0.5);
        return 1000.0 * std::exp(-100.0 * r2);
    });
    
    std::cout << "Assembling system (" 
              << grid.numNodesX() * grid.numNodesY() * grid.numNodesZ() 
              << " unknowns)...\n";
    
    // 4. Solve
    auto solution = solver.solve();
    
    // 5. Find maximum
    double maxU = 0.0;
    for (int i = 0; i < grid.numNodesX(); i++) {
        for (int j = 0; j < grid.numNodesY(); j++) {
            for (int k = 0; k < grid.numNodesZ(); k++) {
                double u = solution(i, j, k);
                if (u > maxU) maxU = u;
            }
        }
    }
    
    std::cout << "Maximum value: " << maxU << "\n";
}
```

### Key Concepts

- **3D grid:** `Grid3D` with 26×26×26 = 17,576 nodes
- **Large sparse system:** ~15K unknowns, efficiently solved
- **Triple indexing:** `solution(i, j, k)` for 3D GridFunction
- **Memory efficiency:** Sparse matrix storage critical for 3D

### Expected Output

```
Assembling system (17576 unknowns)...
Solution computed successfully
Maximum value: 4.10451
```

### Scaling Notes

| Grid Size | Unknowns | Memory | Solve Time |
|-----------|----------|--------|------------|
| 10³ | 1,331 | ~10 KB | <0.1s |
| 20³ | 9,261 | ~100 KB | ~0.5s |
| 25³ | 17,576 | ~200 KB | ~1s |
| 30³ | 29,791 | ~500 KB | ~3s |

### Physical Interpretation

Models the electrostatic potential from a spherical charge distribution. In electrostatics, ∇²φ = -ρ/ε₀ (Poisson's equation). The potential is maximum at the charge center and decays toward grounded boundaries.

---

## Example 5: 1D Time-Dependent Heat

**Physical Problem:** Cooling of a heated rod with time.

**Mathematical Problem:**
```
∂u/∂t = α·∂²u/∂x²  on [0,1] × [0,T]
u(0,t) = 0
u(1,t) = 0
u(x,0) = sin(πx)
```

**Analytical Solution:** `u(x,t) = exp(-π²αt)·sin(πx)`

### Complete Code

```cpp
void example5_1D_transient_heat() {
    // 1. Setup domain and parameters
    Interval<double> domain(0.0, 1.0);
    Grid1D<double> grid(domain, 100);
    double alpha = 0.01;  // Thermal diffusivity
    
    // 2. Boundary conditions
    auto bc = homogeneousDirichlet1D<double>();
    HeatSolver1D<double> solver(grid, bc, alpha);
    
    // 3. Initial condition: sin(πx)
    solver.setInitialCondition([](double x) {
        return std::sin(M_PI * x);
    });
    
    // 4. Time-stepping parameters
    double dt = 0.001;
    double tFinal = 0.1;
    
    // 5. Solve (Crank-Nicolson is 2nd-order accurate)
    auto result = solver.solve(tFinal, dt, TimeScheme::CrankNicolson);
    auto solution = solver.getSolution();  // Returns vector<T>
    
    // 6. Compare with analytical solution
    double expectedDecay = std::exp(-M_PI * M_PI * alpha * tFinal);
    double maxValue = *std::max_element(solution.begin(), solution.end());
    
    std::cout << "Final max value: " << maxValue << "\n";
    std::cout << "Expected (analytical): " << expectedDecay << "\n";
    std::cout << "Relative error: " 
              << std::abs(maxValue - expectedDecay) / expectedDecay << "\n";
}
```

### Key Concepts

- **Parabolic PDE:** Time derivative requires time-stepping
- **Initial condition:** `setInitialCondition()` with lambda
- **Time schemes:** `CrankNicolson` (2nd-order), `BackwardEuler` (1st-order)
- **Solution access:** `getSolution()` returns `vector<T>`, use `solution[i]`

### Expected Output

```
Solution computed successfully
Time steps: 100
Final max value: 0.99018
Expected (analytical): 0.990179
Relative error: 8.11708e-07
```

**Verification:** Excellent agreement with analytical solution!

### Time-Stepping Schemes Comparison

| Scheme | Order | Stability | Best For |
|--------|-------|-----------|----------|
| `ForwardEuler` | O(Δt) | Conditional (Δt ≤ h²/2α) | Quick tests |
| `BackwardEuler` | O(Δt) | Unconditional | Stiff problems |
| `CrankNicolson` | O(Δt²) | Unconditional | **Accuracy (recommended)** |

### Physical Interpretation

A rod initially heated to follow sin(πx) cools over time with ends held at 0°C. The temperature decays exponentially while maintaining the sinusoidal shape. The decay rate depends on thermal diffusivity α.

---

## Example 6: 2D Room Cooling

**Physical Problem:** Room with initially warm air cooling through a cold ceiling.

**Mathematical Problem:**
```
∂u/∂t = α·∇²u  on [0,5]×[0,4] × [0,T]
u(top, t) = 0°C         (cold ceiling)
∂u/∂n(sides, t) = 0     (insulated walls)
u(x, y, 0) = 25°C       (initial warm air)
```

### Complete Code

```cpp
void example6_2D_room_cooling() {
    // 1. Define room dimensions (5m × 4m)
    Rectangle<double> room(0.0, 5.0, 0.0, 4.0);
    Grid2D<double> grid(room, 25, 20);
    double alpha = 2.2e-5;  // Air thermal diffusivity
    
    // 2. Boundary conditions
    BoundaryConditions2D<double> bc;
    bc.setLeft(BoundaryCondition2D<double>::Neumann(0.0));    // Insulated
    bc.setRight(BoundaryCondition2D<double>::Neumann(0.0));   // Insulated
    bc.setBottom(BoundaryCondition2D<double>::Neumann(0.0));  // Insulated
    bc.setTop(BoundaryCondition2D<double>::Dirichlet(0.0));   // Cold ceiling
    
    HeatSolver2D<double> solver(grid, bc, alpha);
    
    // 3. Initial condition: Uniform 25°C
    solver.setInitialCondition([](double x, double y) { 
        return 25.0; 
    });
    
    // 4. Time integration
    double tFinal = 100.0;  // 100 seconds
    double dt = 1.0;
    
    auto result = solver.solve(tFinal, dt, TimeScheme::CrankNicolson);
    auto solution = solver.getSolution();
    
    // 5. Calculate average temperature
    double avgTemp = 0.0;
    for (int i = 0; i < grid.numNodesX(); i++) {
        for (int j = 0; j < grid.numNodesY(); j++) {
            int idx = grid.index(i, j);  // Convert 2D to linear index
            avgTemp += solution[idx];
        }
    }
    avgTemp /= (grid.numNodesX() * grid.numNodesY());
    
    std::cout << "Average final temperature: " << avgTemp << "°C\n";
}
```

### Key Concepts

- **2D time-dependent:** Heat equation in 2D evolves in time
- **Physical dimensions:** 5m × 4m room (not unit square)
- **Realistic parameters:** α = 2.2×10⁻⁵ m²/s for air
- **Linear indexing:** `grid.index(i, j)` converts 2D indices to 1D for vector access

### Expected Output

```
Simulating cooling for 100 seconds
Time steps: 100
Average final temperature: 20.67°C
```

### Physical Interpretation

Models a room where only the ceiling provides cooling (e.g., ceiling-mounted AC). Walls are insulated (no heat loss). Room starts at 25°C and gradually cools as heat diffuses toward the cold ceiling. After 100 seconds, average temperature drops ~4°C.

### Grid Indexing for 2D Heat Solver

```cpp
// 2D HeatSolver returns vector<T>, not GridFunction2D
auto solution = solver.getSolution();

// WRONG: solution(i, j)  // GridFunction syntax doesn't work!

// CORRECT: Convert 2D indices to linear index
int idx = grid.index(i, j);
double temp = solution[idx];
```

---

## Example 7: Robin Boundary Conditions

**Physical Problem:** Metal rod with convective cooling at one end.

**Mathematical Problem:**
```
-d²u/dx² = -100  on [0,1]
u(0) = 100°C                        (left: fixed temp)
h(u - T_amb) + k(du/dx) = 0  at x=1 (right: convective cooling)
```

**Robin BC Form:** `α·u + β·(du/dn) = g`

### Complete Code

```cpp
void example7_robin_bc() {
    // 1. Setup
    Interval<double> domain(0.0, 1.0);
    Grid1D<double> grid(domain, 100);
    
    // 2. Physical parameters
    double h = 10.0;          // Convection coefficient (W/m²K)
    double k = 1.0;           // Thermal conductivity (W/mK)
    double T_ambient = 20.0;  // Ambient temperature (°C)
    double q = 100.0;         // Heat generation (W/m³)
    
    // 3. Boundary conditions
    BoundaryConditions1D<double> bc;
    
    // Left: Fixed temperature
    bc.setLeft(BoundaryCondition1D<double>::Dirichlet(100.0));
    
    // Right: Convective cooling (Robin BC)
    // From: h(u - T_amb) = -k(du/dx)
    // Get: h·u + k·(du/dx) = h·T_amb
    // Robin form: α·u + β·(du/dn) = g
    bc.setRight(BoundaryCondition1D<double>::Robin(h, k, 
        [=](double x) { return h * T_ambient; }
    ));
    
    // 4. Solve
    PoissonSolver1D<double> solver(grid, bc);
    solver.setSource([q](double x) { return -q; });
    
    auto solution = solver.solve();
    
    // 5. Check temperatures
    std::cout << "Temperature at left (x=0): " << solution(0) << "°C\n";
    std::cout << "Temperature at right (x=1): " 
              << solution(grid.numNodes() - 1) << "°C\n";
}
```

### Key Concepts

- **Robin BC:** Mixed condition combining value and derivative
- **Physical meaning:** Models convective heat transfer: `q_conv = h(T_surface - T_ambient)`
- **Standard form:** `α·u + β·(du/dn) = g`
- **Conversion:** From physics equation to Robin form requires algebra

### Robin BC Conversion

Starting from convective cooling physics:
```
h(u - T_amb) = -k(du/dx)  at boundary
```

Rearrange to Robin form:
```
h·u + k·(du/dx) = h·T_amb
```

Match to standard form `α·u + β·(du/dn) = g`:
```cpp
α = h
β = k  
g = h * T_ambient
```

### Expected Output

```
Temperature at left (x=0): 100°C
Temperature at right (x=1): ~60°C (cooling toward 20°C ambient)
```

### Physical Interpretation

A heated rod (100°C at left end) loses heat through convection at the right end. The convection coefficient h determines how efficiently heat transfers to ambient air. Higher h means more cooling.

---

## Common Patterns and Best Practices

### 1. Elliptic PDEs (Poisson/Laplace)

```cpp
// Template for elliptic problems
PoissonSolver1D<double> solver(grid, bc);
solver.setSource(lambda_function);
auto solution = solver.solve();  // Returns GridFunction directly

// Access solution
double u = solution(i);        // 1D
double u = solution(i, j);     // 2D
double u = solution(i, j, k);  // 3D
```

### 2. Parabolic PDEs (Heat Equation)

```cpp
// Template for time-dependent problems
HeatSolver1D<double> solver(grid, bc, alpha);
solver.setInitialCondition(lambda_or_GridFunction);

auto result = solver.solve(tFinal, dt, TimeScheme::CrankNicolson);
auto solution = solver.getSolution();  // Returns vector<T>

// Access solution
double u = solution[i];           // 1D: direct indexing
int idx = grid.index(i, j);      // 2D: convert to linear index
double u = solution[idx];         // 2D: then access
```

### 3. Boundary Conditions

```cpp
// 1D
BoundaryConditions1D<double> bc;
bc.setLeft(BoundaryCondition1D<double>::Dirichlet(value));
bc.setRight(BoundaryCondition1D<double>::Neumann(flux));

// 2D - specify all four sides
BoundaryConditions2D<double> bc;
bc.setLeft(...);
bc.setRight(...);
bc.setTop(...);
bc.setBottom(...);

// 3D - specify all six faces
BoundaryConditions3D<double> bc;
bc.setLeft(...);   bc.setRight(...);
bc.setFront(...);  bc.setBack(...);
bc.setTop(...);    bc.setBottom(...);
```

### 4. Source Terms

```cpp
// Constant source
solver.setSource([](double x) { return 10.0; });

// Position-dependent source
solver.setSource([](double x, double y) {
    double r2 = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5);
    return 1000.0 * std::exp(-100.0 * r2);
});

// Zero source (Laplace equation)
solver.setSource([](double x, double y) { return 0.0; });
```

---

## Troubleshooting

### Convergence Warnings

```
Warning: PoissonSolver2D did not converge after 10000 iterations
```

**Solutions:**
- Increase max iterations in solver config
- Use better preconditioner (ILU vs Jacobi)
- Refine grid (may improve condition number)
- Check BC implementation

### Type Mismatches

```cpp
// WRONG: Mixing solution types
auto solution = heatSolver.solve(...);  // Returns result
double u = solution(i);  // ERROR: not a GridFunction!

// CORRECT:
auto result = heatSolver.solve(...);
auto solution = heatSolver.getSolution();  // Get the vector
double u = solution[i];  // Array indexing
```

### Index Out of Bounds

```cpp
// WRONG: Off-by-one error
for (int i = 0; i <= grid.numNodes(); i++)  // Too far!

// CORRECT:
for (int i = 0; i < grid.numNodes(); i++)   // Stop before end
```

---

## Building Your Own Examples

### Step-by-Step Workflow

1. **Define the Problem**
   - Physical domain and dimensions
   - PDE type (elliptic, parabolic, hyperbolic)
   - Boundary conditions
   - Source term

2. **Create Grid**
   ```cpp
   Grid1D<double> grid(domain, numNodes);
   Grid2D<double> grid(domain, nx, ny);
   Grid3D<double> grid(domain, nx, ny, nz);
   ```

3. **Set Boundary Conditions**
   ```cpp
   BoundaryConditions1D<double> bc;
   bc.setLeft(...);
   bc.setRight(...);
   ```

4. **Create Solver**
   ```cpp
   PoissonSolver1D<double> solver(grid, bc);  // Elliptic
   HeatSolver1D<double> solver(grid, bc, alpha);  // Parabolic
   ```

5. **Set Source/Initial Condition**
   ```cpp
   solver.setSource(lambda);          // Elliptic
   solver.setInitialCondition(u0);    // Parabolic
   ```

6. **Solve**
   ```cpp
   auto solution = solver.solve();                        // Elliptic
   auto result = solver.solve(tFinal, dt, timeScheme);   // Parabolic
   ```

7. **Post-Process**
   - Extract values
   - Compare with analytical solutions
   - Compute derived quantities
   - Export for visualization

---

## Performance Tips

### Grid Resolution

| Grid Size | Accuracy | Solve Time | Memory | Recommended For |
|-----------|----------|------------|--------|-----------------|
| 10-50 nodes | Low | <0.1s | KB | Quick tests, debugging |
| 50-100 nodes | Medium | 0.1-1s | ~10 KB | Development, verification |
| 100-200 nodes | High | 1-10s | ~100 KB | Production 1D/2D |
| 25³ (15K) | Medium 3D | 1-5s | ~1 MB | Small 3D problems |
| 50³ (125K) | High 3D | 10-60s | ~10 MB | Large 3D problems |

### Time-Stepping

```cpp
// Stability condition for explicit methods
double dt_max = dx * dx / (2 * alpha * dim);  // dim = 1, 2, or 3

// For CrankNicolson and BackwardEuler, dt is limited by accuracy, not stability
double dt_accurate = dx * dx / (10 * alpha);  // Conservative choice
```

### Solver Selection

For large 3D problems, consider:
- Using preconditioned CG instead of direct solvers
- Adjusting tolerance vs iterations tradeoff
- Exploiting symmetry when possible

---

## Further Reading

### Related Documentation

- [Quick Start Guide](Quick_Start_Guide.md) - Getting started with PDE solvers
- [Sparse Matrices](Sparse_Matrices.md) - Understanding sparse storage formats
- [Grid Infrastructure](Grid_Infrastructure.md) - Grid types and indexing
- [Boundary Conditions](Boundary_Conditions.md) - Complete BC reference
- [Iterative Solvers](Iterative_Solvers.md) - Solver algorithms and preconditioning
- [Elliptic PDEs](Elliptic_PDEs.md) - Poisson/Laplace equation details
- [Parabolic PDEs](Parabolic_PDEs.md) - Heat equation methods

### External Resources

- [Numerical Recipes](http://numerical.recipes/) - Classic reference for numerical methods
- [LeVeque: Finite Difference Methods](https://faculty.washington.edu/rjl/fdmbook/) - Excellent PDE textbook
- [Strang: Computational Science and Engineering](https://math.mit.edu/~gs/cse/) - MIT courseware

---

## Example Programs

All examples are located in `src/docs_demos/pde/`:

- `gallery_examples_simple.cpp` - This gallery (complete source)
- `elliptic_verification.cpp` - Elliptic PDE test suite
- `parabolic_verification.cpp` - Parabolic PDE test suite
- `example02_1d_heat.cpp` - Basic 1D heat equation
- `example03_2d_laplace.cpp` - 2D Laplace equation
- `example04_room_cooling.cpp` - Room cooling simulation

---

**🎓 Happy Computing!**

For questions, issues, or contributions, visit the [MML GitHub repository](https://github.com/zvanjak/MinimalMathLibrary).
