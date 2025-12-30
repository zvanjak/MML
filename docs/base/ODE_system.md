# ODE Systems

**Sources:**  
- `mml/interfaces/IODESystem.h` - ODE system interfaces
- `mml/base/ODESystem.h` - ODE system classes
- `mml/base/ODESystemSolution.h` - Solution storage and post-processing

Classes for representing and solving systems of ordinary differential equations (ODEs). MML provides a clean separation between defining ODE systems and solving them, allowing flexible solver choice and solution post-processing.

## Overview

**Core Concepts:**
- **ODE System**: Defines dynamics **x'** = **f**(t, **x**) for state vector **x**
- **Initial Conditions**: Starting state **x**(t₀) = **x**₀
- **Solution**: Time evolution **x**(t) stored as discrete steps
- **Post-processing**: Convert solutions to interpolated functions or parametric curves

**System Types:**
- **Basic ODE System**: Provides derivative function **f**(t, **x**)
- **System with Jacobian**: Includes analytical Jacobian ∂**f**/∂**x** for stiff problems
- **Parametrized System**: Systems with adjustable parameters

**Solution Storage:**
- Automatic resizing as integration proceeds
- Efficient storage with dynamic extension
- Export to interpolated functions or parametric curves
- Statistics tracking (successful/failed steps)

---

## Quick Reference

### ODE System Classes

| Class | Purpose | Use Case |
|-------|---------|----------|
| `IODESystem` | Interface | Define custom ODE systems |
| `ODESystem` | Basic implementation | Non-stiff systems |
| `IODESystemWithJacobian` | Interface with Jacobian | Stiff systems (interface) |
| `ODESystemWithJacobian` | Implementation with Jacobian | Stiff systems with analytical Jacobian |
| `IODESystemParametrized` | Parametrized interface | Systems with adjustable parameters |

### Solution Class

| Class | Purpose |
|-------|---------|
| `ODESystemSolution` | Stores solution trajectory, provides interpolation/export |

---

## ODE System Interfaces

### IODESystem

Base interface for all ODE systems.

**Interface:**
```cpp
class IODESystem
{
public:
    virtual int getDim() const = 0;
    virtual void derivs(const Real t, const Vector<Real>& x, 
                        Vector<Real>& dxdt) const = 0;
    
    virtual std::string getVarName(int ind) const;
};
```

**Methods:**
- `getDim()`: Return system dimension (number of equations)
- `derivs(t, x, dxdt)`: Compute derivatives dxdt = f(t, x)
- `getVarName(ind)`: Optional variable names for output/visualization

**Mathematical Definition:**
```
System: x' = f(t, x)
where:
  x ∈ ℝⁿ  (state vector)
  t ∈ ℝ   (time)
  f: ℝ × ℝⁿ → ℝⁿ  (derivative function)
```

### IODESystemWithJacobian

Interface for systems that can provide analytical Jacobian.

**Interface:**
```cpp
class IODESystemWithJacobian : public IODESystem
{
public:
    virtual void jacobian(const Real t, const Vector<Real>& x, 
                          Vector<Real>& dxdt, 
                          Matrix<Real>& dydx) const = 0;
};
```

**Jacobian Matrix:**
```
J = ∂f/∂x = [∂fᵢ/∂xⱼ]

For stiff systems, providing analytical Jacobian improves:
- Stability of implicit methods
- Step size control
- Convergence of Newton iterations
```

**When to Use:**
- Stiff differential equations
- Implicit integration methods
- When analytical derivatives are known

### IODESystemParametrized

Interface for systems with adjustable parameters.

**Interface:**
```cpp
class IODESystemParametrized : public IODESystem
{
public:
    virtual int getNumParam() const = 0;
    virtual Real getParam(int i) const = 0;
    virtual void setParam(int i, Real val) = 0;
    
    virtual Vector<Real> getParams() const = 0;
    virtual void setParams(const Vector<Real>&) = 0;
};
```

**Use Cases:**
- Parameter studies
- Optimization problems
- Sensitivity analysis
- System identification

---

## ODE System Implementations

### ODESystem

Basic ODE system implementation using function pointer.

**Constructor:**
```cpp
ODESystem();  // Empty system
ODESystem(int n, void (*func)(Real, const Vector<Real>&, Vector<Real>&));
```

**Parameters:**
- `n`: System dimension
- `func`: Derivative function with signature `void func(Real t, const Vector<Real>& x, Vector<Real>& dxdt)`

**Example 1: Exponential Decay**
```cpp
// System: x' = -k*x
void exponential_decay(Real t, const Vector<Real>& x, Vector<Real>& dxdt)
{
    Real k = 0.5;  // Decay rate
    dxdt[0] = -k * x[0];
}

ODESystem decay_system(1, exponential_decay);

// Solve with initial condition x(0) = 1.0
Vector<Real> x0(1);
x0[0] = 1.0;

// Use with solver (see algorithms documentation)
```

**Example 2: Harmonic Oscillator**
```cpp
// System: x'' + ω²x = 0
// Rewrite as: x₁' = x₂, x₂' = -ω²x₁
void harmonic_oscillator(Real t, const Vector<Real>& x, Vector<Real>& dxdt)
{
    Real omega = 2.0 * Constants::PI;  // Frequency
    
    dxdt[0] = x[1];              // x' = v
    dxdt[1] = -omega*omega*x[0]; // v' = -ω²x
}

ODESystem oscillator(2, harmonic_oscillator);

// Initial conditions: x(0) = 1.0, v(0) = 0.0
Vector<Real> x0{1.0, 0.0};
```

**Example 3: Lorenz System**
```cpp
// Chaotic system: 
//   x' = σ(y - x)
//   y' = x(ρ - z) - y
//   z' = xy - βz
void lorenz_system(Real t, const Vector<Real>& x, Vector<Real>& dxdt)
{
    Real sigma = 10.0;
    Real rho = 28.0;
    Real beta = 8.0/3.0;
    
    dxdt[0] = sigma * (x[1] - x[0]);
    dxdt[1] = x[0] * (rho - x[2]) - x[1];
    dxdt[2] = x[0] * x[1] - beta * x[2];
}

ODESystem lorenz(3, lorenz_system);

// Initial conditions
Vector<Real> x0{1.0, 1.0, 1.0};
```

### ODESystemWithJacobian

ODE system with analytical Jacobian for stiff problems.

**Constructor:**
```cpp
ODESystemWithJacobian();  // Empty
ODESystemWithJacobian(int n,
    void (*func)(Real, const Vector<Real>&, Vector<Real>&),
    void (*jac)(const Real, const Vector<Real>&, Vector<Real>&, Matrix<Real>&)
);
```

**Parameters:**
- `n`: System dimension
- `func`: Derivative function
- `jac`: Jacobian function with signature `void jac(Real t, const Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& J)`

**Example: Stiff Chemical Reaction**
```cpp
// Robertson problem (stiff):
//   y₁' = -0.04*y₁ + 10⁴*y₂*y₃
//   y₂' = 0.04*y₁ - 10⁴*y₂*y₃ - 3×10⁷*y₂²
//   y₃' = 3×10⁷*y₂²

void robertson_derivs(Real t, const Vector<Real>& y, Vector<Real>& dydt)
{
    dydt[0] = -0.04 * y[0] + 1.0e4 * y[1] * y[2];
    dydt[1] = 0.04 * y[0] - 1.0e4 * y[1] * y[2] - 3.0e7 * y[1] * y[1];
    dydt[2] = 3.0e7 * y[1] * y[1];
}

void robertson_jacobian(const Real t, const Vector<Real>& y, 
                        Vector<Real>& dydt, Matrix<Real>& J)
{
    // J[i][j] = ∂fᵢ/∂yⱼ
    J[0][0] = -0.04;
    J[0][1] = 1.0e4 * y[2];
    J[0][2] = 1.0e4 * y[1];
    
    J[1][0] = 0.04;
    J[1][1] = -1.0e4 * y[2] - 6.0e7 * y[1];
    J[1][2] = -1.0e4 * y[1];
    
    J[2][0] = 0.0;
    J[2][1] = 6.0e7 * y[1];
    J[2][2] = 0.0;
}

ODESystemWithJacobian robertson(3, robertson_derivs, robertson_jacobian);

// Initial conditions
Vector<Real> y0{1.0, 0.0, 0.0};

// Requires stiff solver (e.g., implicit Runge-Kutta, BDF)
```

---

## ODE System Solution

### ODESystemSolution

Stores and post-processes ODE solution trajectory.

**Constructors:**
```cpp
ODESystemSolution(Real t1, Real t2, int dim, int maxSteps);
ODESystemSolution(Real t1, Real t2, int dim);  // Default maxSteps=1000
```

**Parameters:**
- `t1`, `t2`: Time interval [t₁, t₂]
- `dim`: System dimension
- `maxSteps`: Initial storage capacity (auto-extends if needed)

### Storage Management

**Internal Storage:**
- `_tval`: Vector of time values
- `_xval`: Matrix of state values (dim × numSteps)
- Automatic resizing when capacity exceeded

**Methods:**
```cpp
void fillValues(int ind, Real t, Vector<Real>& x);  // Store solution step
void setFinalSize(int numSteps);                     // Trim to actual size
```

### Solution Access

**Retrieve Values:**
```cpp
Vector<Real> getTValues() const;              // All time values
Matrix<Real> getXValues() const;              // All state values
Vector<Real> getXValues(int component) const; // Single component
Vector<Real> getXValuesAtEnd() const;         // Final state
```

**Statistics:**
```cpp
int getNumStepsOK() const;       // Successful steps
int getNumStepsBad() const;      // Rejected steps
int getTotalNumSteps() const;    // Total attempted
int getTotalSavedSteps() const;  // Actually stored
```

**Example:**
```cpp
ODESystemSolution solution(0.0, 10.0, 3);

// After solving...
Vector<Real> times = solution.getTValues();
Vector<Real> x_component = solution.getXValues(0);  // First component

std::cout << "Successful steps: " << solution.getNumStepsOK() << "\n";
std::cout << "Failed steps: " << solution.getNumStepsBad() << "\n";
```

### Solution Export

#### Interpolated Functions

Export individual components as interpolated functions.

**Methods:**
```cpp
LinearInterpRealFunc getSolAsLinInterp(int component) const;
PolynomInterpRealFunc getSolAsPolyInterp(int component, int polyOrder) const;
SplineInterpRealFunc getSolAsSplineInterp(int component) const;
```

**Example:**
```cpp
// Solve ODE system
ODESystemSolution solution = solver.solve(system, x0, 0.0, 10.0);

// Export first component as spline
SplineInterpRealFunc x_interp = solution.getSolAsSplineInterp(0);

// Evaluate at arbitrary times
Real x_at_5 = x_interp(5.0);
Real x_at_7_5 = x_interp(7.5);

// Use with other MML functions
Real deriv = Derivation::Derive(x_interp, 5.0, nullptr);
```

#### Parametric Curves

Convert phase space trajectories to parametric curves.

**Methods:**
```cpp
SplineInterpParametricCurve<2> getSolAsParamCurve2D(int ind1, int ind2) const;
SplineInterpParametricCurve<3> getSolAsParamCurve3D(int ind1, int ind2, int ind3) const;
```

**Example: Phase Space Trajectory**
```cpp
// Solve oscillator
ODESystemSolution solution = solver.solve(oscillator, x0, 0.0, 20.0);

// Extract phase space curve (position vs velocity)
auto phase_curve = solution.getSolAsParamCurve2D(0, 1);  // x vs v

// Analyze curve
Real arc_length = PathIntegration::ParametricCurveLength(phase_curve, 0.0, 20.0);

// Visualize
Visualizer::VisualizeParamCurve2D(phase_curve, "phase_space", 0.0, 20.0, 500);
```

**Example: 3D Trajectory (Lorenz Attractor)**
```cpp
// Solve Lorenz system
ODESystemSolution solution = solver.solve(lorenz, x0, 0.0, 50.0);

// Extract 3D trajectory
auto trajectory = solution.getSolAsParamCurve3D(0, 1, 2);

// Visualize attractor
Visualizer::VisualizeParamCurve3D(trajectory, "lorenz_attractor", 0.0, 50.0, 1000);
```

---

## Complete Examples

### Example 1: Simple Pendulum

```cpp
#include "base/ODESystem.h"
#include "algorithms/ODESolver.h"

void Example_Simple_Pendulum()
{
    // Pendulum: θ'' + (g/L)sin(θ) = 0
    // System: θ' = ω, ω' = -(g/L)sin(θ)
    
    auto pendulum = [](Real t, const Vector<Real>& x, Vector<Real>& dxdt) {
        Real g = 9.81;   // Gravity
        Real L = 1.0;    // Length
        
        dxdt[0] = x[1];                    // θ' = ω
        dxdt[1] = -(g/L) * std::sin(x[0]); // ω' = -(g/L)sin(θ)
    };
    
    ODESystem system(2, pendulum);
    
    // Initial conditions: θ=π/4 (45°), ω=0
    Vector<Real> x0{Constants::PI/4, 0.0};
    
    // Solve using RK4
    ODESystemSolver<RK4Step> solver(system);
    auto solution = solver.solve(0.0, 10.0, x0, 0.01);
    
    // Extract angle as function of time
    auto theta_t = solution.getSolAsSplineInterp(0);
    
    // Phase space trajectory
    auto phase_curve = solution.getSolAsParamCurve2D(0, 1);
    
    // Verify energy conservation
    Real E0 = 0.5 * x0[1]*x0[1] + (9.81/1.0) * (1 - std::cos(x0[0]));
    
    Vector<Real> x_final = solution.getXValuesAtEnd();
    Real E_final = 0.5 * x_final[1]*x_final[1] + (9.81/1.0) * (1 - std::cos(x_final[0]));
    
    std::cout << "Energy drift: " << std::abs(E_final - E0) << "\n";
}
```

### Example 2: N-Body Problem

```cpp
void Example_TwoBody_System()
{
    // Two-body gravitational system
    // State: [x₁, y₁, vx₁, vy₁, x₂, y₂, vx₂, vy₂]
    
    auto two_body = [](Real t, const Vector<Real>& x, Vector<Real>& dxdt) {
        Real G = 6.67430e-11;  // Gravitational constant
        Real m1 = 1.0e24;       // Mass 1
        Real m2 = 1.0e24;       // Mass 2
        
        // Positions
        Real dx = x[4] - x[0];  // x₂ - x₁
        Real dy = x[5] - x[1];  // y₂ - y₁
        Real r = std::sqrt(dx*dx + dy*dy);
        Real r3 = r*r*r;
        
        // Forces
        Real Fx = G * m1 * m2 * dx / r3;
        Real Fy = G * m1 * m2 * dy / r3;
        
        // Derivatives
        dxdt[0] = x[2];           // vx₁
        dxdt[1] = x[3];           // vy₁
        dxdt[2] = Fx / m1;        // ax₁
        dxdt[3] = Fy / m1;        // ay₁
        
        dxdt[4] = x[6];           // vx₂
        dxdt[5] = x[7];           // vy₂
        dxdt[6] = -Fx / m2;       // ax₂
        dxdt[7] = -Fy / m2;       // ay₂
    };
    
    ODESystem system(8, two_body);
    
    // Initial conditions: circular orbit
    Vector<Real> x0(8);
    x0[0] = 1.0e8;  x0[1] = 0.0;      // Body 1 position
    x0[2] = 0.0;    x0[3] = 1000.0;   // Body 1 velocity
    x0[4] = -1.0e8; x0[5] = 0.0;      // Body 2 position
    x0[6] = 0.0;    x0[7] = -1000.0;  // Body 2 velocity
    
    // Solve
    ODESystemSolver<RK4Step> solver(system);
    auto solution = solver.solve(0.0, 1000.0, x0, 1.0);
    
    // Extract trajectories
    auto body1_traj = solution.getSolAsParamCurve2D(0, 1);
    auto body2_traj = solution.getSolAsParamCurve2D(4, 5);
}
```

### Example 3: Stiff Van der Pol Oscillator

```cpp
void Example_Stiff_VanDerPol()
{
    // Van der Pol: x'' - μ(1-x²)x' + x = 0
    // For large μ, system is stiff
    
    Real mu = 1000.0;  // Stiffness parameter
    
    auto van_der_pol_derivs = [mu](Real t, const Vector<Real>& x, Vector<Real>& dxdt) {
        dxdt[0] = x[1];
        dxdt[1] = mu * (1.0 - x[0]*x[0]) * x[1] - x[0];
    };
    
    auto van_der_pol_jacobian = [mu](const Real t, const Vector<Real>& x, 
                                      Vector<Real>& dxdt, Matrix<Real>& J) {
        J[0][0] = 0.0;
        J[0][1] = 1.0;
        J[1][0] = -2.0 * mu * x[0] * x[1] - 1.0;
        J[1][1] = mu * (1.0 - x[0]*x[0]);
    };
    
    ODESystemWithJacobian system(2, van_der_pol_derivs, van_der_pol_jacobian);
    
    Vector<Real> x0{2.0, 0.0};
    
    // Requires stiff solver (implicit method)
    // StiffSolver solver(system);
    // auto solution = solver.solve(0.0, 3000.0, x0);
    
    // Phase portrait
    // auto phase_curve = solution.getSolAsParamCurve2D(0, 1);
}
```

### Example 4: Parametrized System Study

```cpp
class ParametrizedOscillator : public IODESystemParametrized
{
private:
    Real _omega;  // Frequency parameter
    Real _damping; // Damping parameter
    
public:
    ParametrizedOscillator(Real omega, Real damping) 
        : _omega(omega), _damping(damping) {}
    
    int getDim() const override { return 2; }
    
    void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
    {
        dxdt[0] = x[1];
        dxdt[1] = -_omega*_omega*x[0] - 2.0*_damping*x[1];
    }
    
    int getNumParam() const override { return 2; }
    
    Real getParam(int i) const override {
        return (i == 0) ? _omega : _damping;
    }
    
    void setParam(int i, Real val) override {
        if (i == 0) _omega = val;
        else _damping = val;
    }
    
    Vector<Real> getParams() const override {
        return Vector<Real>{_omega, _damping};
    }
    
    void setParams(const Vector<Real>& params) override {
        _omega = params[0];
        _damping = params[1];
    }
};

void Example_Parameter_Study()
{
    Vector<Real> x0{1.0, 0.0};
    
    // Study damping effect
    for (Real damping = 0.0; damping <= 2.0; damping += 0.2)
    {
        ParametrizedOscillator system(2.0, damping);
        
        ODESystemSolver<RK4Step> solver(system);
        auto solution = solver.solve(0.0, 10.0, x0, 0.01);
        
        // Analyze decay rate
        auto x_t = solution.getSolAsSplineInterp(0);
        
        std::cout << "Damping: " << damping 
                  << ", Final amplitude: " << std::abs(x_t(10.0)) << "\n";
    }
}
```

---

## Best Practices

### System Definition

**1. Choose appropriate representation:**
```cpp
// Simple: Function pointer
ODESystem simple(n, derivs_func);

// Stiff: Include Jacobian
ODESystemWithJacobian stiff(n, derivs_func, jacobian_func);

// Complex: Derive from interface
class CustomSystem : public IODESystem { /* ... */ };
```

**2. Validate dimensions:**
```cpp
void safe_derivs(Real t, const Vector<Real>& x, Vector<Real>& dxdt)
{
    assert(x.size() == expected_dim);
    assert(dxdt.size() == expected_dim);
    // Compute derivatives...
}
```

**3. Document physics:**
```cpp
// Good: Clear physical meaning
void damped_oscillator(Real t, const Vector<Real>& x, Vector<Real>& dxdt)
{
    // x[0] = position, x[1] = velocity
    // Parameters: ω = natural frequency, γ = damping
    Real omega = 2.0;
    Real gamma = 0.1;
    
    dxdt[0] = x[1];                           // v
    dxdt[1] = -omega*omega*x[0] - 2*gamma*x[1];  // a
}
```

### Solution Storage

**1. Estimate storage needs:**
```cpp
// Conservative estimate
int max_steps = static_cast<int>((t2 - t1) / min_step_size) + 100;
ODESystemSolution solution(t1, t2, dim, max_steps);
```

**2. Trim after solving:**
```cpp
auto solution = solver.solve(/* ... */);
solution.setFinalSize(actual_steps);  // Reduce memory
```

### Post-Processing

**1. Check solution quality:**
```cpp
int failed = solution.getNumStepsBad();
int total = solution.getTotalNumSteps();

if (failed > 0.1 * total)
    std::cerr << "Warning: " << failed << " rejected steps!\n";
```

**2. Use appropriate interpolation:**
```cpp
// Smooth solutions: Spline
auto smooth = solution.getSolAsSplineInterp(component);

// Discontinuous: Linear
auto discontinuous = solution.getSolAsLinInterp(component);
```

---

## Integration with MML

### With ODE Solvers (algorithms/)

```cpp
// Fixed-step solvers
ODESystemSolver<RK4Step> rk4_solver(system);
ODESystemSolver<DormandPrince54> dp54_solver(system);

// Adaptive solvers
AdaptiveODESolver adaptive(system, tolerance);

// Specialized solvers
LeapfrogSolver symplectic(hamiltonian_system);
```

### With Interpolation

```cpp
// Solution → interpolated function
auto x_interp = solution.getSolAsSplineInterp(0);

// Use in further calculations
Real derivative = Derivation::Derive(x_interp, t, nullptr);
Real integral = Integration::Integrate(x_interp, t1, t2);
```

### With Curves and Surfaces

```cpp
// 2D phase portrait
auto phase_curve = solution.getSolAsParamCurve2D(0, 1);
Real curvature = phase_curve.getCurvature(t);

// 3D trajectory
auto trajectory = solution.getSolAsParamCurve3D(0, 1, 2);
Vec3Cart tangent = trajectory.getTangent(t);
```

### With Visualization

```cpp
// Plot time series
Vector<Real> times = solution.getTValues();
Vector<Real> values = solution.getXValues(0);
Visualizer::PlotData(times, values, "solution");

// Plot phase space
auto curve = solution.getSolAsParamCurve2D(0, 1);
Visualizer::VisualizeParamCurve2D(curve, "phase", t1, t2, 500);
```

---

## Runnable Examples

| Example | Source File | Description |
|---------|-------------|-------------|
| ODE Solvers Demo | [docs_demo_ode_solvers.cpp](../../src/docs_demos/docs_demo_ode_solvers.cpp) | Comprehensive ODE solver demonstrations |

### Demo Functions

**IODESystem Interface:**
- `Docs_Demo_IODESystem_Interface()` - Interface methods, predefined systems
- `Docs_Demo_ODESystem_Creation()` - Creating ODE systems, TestBeds examples

**Step Calculators & Solvers:**
- `Docs_Demo_Step_Calculators()` - Euler, RK4, Leapfrog comparison on SHO
- `Docs_Demo_Adaptive_Integrators()` - Cash-Karp, Dormand-Prince integrators
- `Docs_Demo_Legacy_Solver_Interface()` - Backward-compatible solver API

**ODESystemSolution Post-Processing:**
- `Docs_Demo_ODESystemSolution_Basics()` - Accessing solution data, statistics
- `Docs_Demo_ODESystemSolution_Interpolation()` - Linear, polynomial, spline interpolation
- `Docs_Demo_ODESystemSolution_ParametricCurve()` - 2D/3D phase space curves

**Physics Examples:**
- `Docs_Demo_Physics_SHO()` - Simple harmonic oscillator with analytical verification
- `Docs_Demo_Physics_Damped_Oscillator()` - Damped oscillator amplitude decay
- `Docs_Demo_Physics_Lorenz()` - Lorenz attractor, sensitive dependence demo
- `Docs_Demo_Physics_Predator_Prey()` - Lotka-Volterra predator-prey dynamics

**Original Demo:**
- `Docs_Demo_ODE_solvers_RungeKutta_4th_order()` - Fixed and adaptive RK methods


---

## See Also

- [Functions.md](Functions.md) - Function interfaces
- [Interpolated_functions.md](Interpolated_functions.md) - Solution interpolation
- [Curves_and_surfaces.md](Curves_and_surfaces.md) - Parametric curves
- [../algorithms/ODE_solvers.md](../algorithms/ODE_solvers.md) - Integration methods
- [Derivation.md](Derivation.md) - Numerical Jacobians
- [../tools/Visualizers_README.md](../tools/Visualizers_README.md) - Solution visualization

---

**References:**
- Hairer, Nørsett, Wanner, *Solving Ordinary Differential Equations I* (2nd ed.), Springer, 1993
- Press et al., *Numerical Recipes* (3rd ed.), Cambridge, 2007
- Ascher, Petzold, *Computer Methods for Ordinary Differential Equations*, SIAM, 1998
