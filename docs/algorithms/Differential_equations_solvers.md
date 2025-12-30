# Differential Equations Solvers

Comprehensive toolkit for solving **Ordinary Differential Equations (ODEs)** - from simple Euler to sophisticated adaptive Runge-Kutta methods.

## Overview

**Solves**: First-order ODE systems of form **dy/dt = f(t, y)** where y ∈ ℝⁿ

**Architecture**:
- **Step Calculators**: Compute single step y(t+h) from y(t) - the numerical "engine"
- **Steppers**: Manage adaptive timestep control for accuracy/efficiency
- **Solvers**: Orchestrate integration from t₁ to t₂, collect results

**Design Philosophy**: Separation of concerns - swap calculators/steppers independently for different problems.

## Quick Reference

### Fixed-Step Methods (ODESystemFixedStepSolver)

Simple, predictable, fixed timestep h = (t₂-t₁)/numSteps

| Calculator | Order | Best For | Pros | Cons |
|-----------|-------|----------|------|------|
| **EulerStep** | 1 | Teaching, quick tests | Simple, fast | Very inaccurate, unstable |
| **Midpoint** | 2 | Low accuracy OK | Better than Euler | Still poor |
| **RungeKutta4 (RK4)** | 4 | Standard workhorse | Good balance | Fixed step only |
| **VelocityVerlet** | 2 | Hamiltonian systems | Energy-conserving | Position/velocity split required |
| **Leapfrog** | 2 | N-body, celestial mechanics | Symplectic, stable | Position/velocity split required |

**Use when**: Simple problems, predictable behavior, benchmarking, teaching

### Adaptive-Step Methods (ODESystemSolver<Stepper>)

Automatic timestep adjustment for accuracy vs efficiency

| Stepper | Order | Error Control | Best For | Performance |
|---------|-------|---------------|----------|-------------|
| **RK5_CashKarp** | 5 | Embedded 4th-order | General non-stiff | ⭐⭐⭐⭐ |
| **DormandPrince5** | 5 | Embedded 4th-order | Production code (non-stiff) | ⭐⭐⭐⭐⭐ |
| **DormandPrince8** | 8 | Embedded 7th-order | High precision needed | ⭐⭐⭐⭐ (slower) |

**Use when**: Accuracy critical, unknown behavior, production code, efficiency matters

### Special-Purpose

| Class | Purpose | Use Case |
|-------|---------|----------|
| **ODESystemLeapfrogSolver** | Hamiltonian dynamics | N-body simulations, celestial mechanics |

---

## Problem Classification

### Stiffness

**Non-Stiff**: Solution components vary on similar time scales
- **Example**: Simple harmonic oscillator, projectile motion
- **Use**: RK4, Dormand-Prince 5/8
- **Characteristic**: All eigenvalues of Jacobian have similar magnitude

**Stiff**: Solution components vary on vastly different time scales
- **Example**: Chemical reactions, electrical circuits with fast/slow dynamics
- **Use**: Implicit methods (not yet implemented), small timesteps with RK4
- **Characteristic**: Jacobian eigenvalues span many orders of magnitude
- **Problem**: Explicit methods need tiny h to maintain stability

### System Structure

**General**: dy/dt = f(t, y)
- Use any method

**Hamiltonian**: H(p,q) conserved, symplectic structure
- **Use**: Leapfrog, Velocity Verlet (preserve energy, symplectic form)
- **Example**: Planetary orbits, molecular dynamics
- **Structure**: y = [q, p], dq/dt = ∂H/∂p, dp/dt = -∂H/∂q

**Second-Order**: d²y/dt² = f(t, y, dy/dt)
- **Convert to first-order**: y₁ = y, y₂ = dy/dt → dy₁/dt = y₂, dy₂/dt = f(t,y₁,y₂)
- Then use any first-order solver

---

## Mathematical Background

### Euler Method

**Simplest**: First-order Taylor expansion

**Formula**:
```
y_{n+1} = y_n + h·f(t_n, y_n)
```

**Error**: Local error ~ O(h²), global error ~ O(h)

**Stability**: Unstable for many problems (especially stiff), needs very small h

### Runge-Kutta Methods

**Idea**: Use **multiple evaluations** within interval to achieve higher order

**RK4** (Classic 4th-order):
```
k₁ = f(t_n, y_n)
k₂ = f(t_n + h/2, y_n + h·k₁/2)
k₃ = f(t_n + h/2, y_n + h·k₂/2)
k₄ = f(t_n + h, y_n + h·k₃)

y_{n+1} = y_n + (h/6)(k₁ + 2k₂ + 2k₃ + k₄)
```

**Error**: Local ~ O(h⁵), global ~ O(h⁴)

**Evaluations**: 4 per step

### Embedded RK Methods (Adaptive)

**Idea**: Compute **two estimates** (orders p and p-1) from same k evaluations
- **Higher-order** estimate used for solution
- **Lower-order** estimate used to estimate error
- **Adjust h** based on error estimate

**Cash-Karp RK5**:
- 5th-order solution
- 4th-order error estimate
- 6 function evaluations per step

**Dormand-Prince RK5(4)**:
- 5th-order solution
- 4th-order embedded error
- 7 function evaluations
- **FSAL property** (First Same As Last): k₁ of next step = k₇ of current step
- Most popular general-purpose ODE solver

**Dormand-Prince RK8(7)**:
- 8th-order solution
- 7th-order embedded error
- 13 function evaluations
- For high-precision requirements

### Step Size Control

**Goal**: Maintain error per step ~ ε (tolerance)

**Algorithm**:
```
1. Estimate local error: err = ||y_{high} - y_{low}||
2. Compute error ratio: r = ε / err
3. If r ≥ 1: Accept step
4. If r < 1: Reject step, reduce h
5. Compute next h:
   h_new = SAFETY × h × r^(1/ORDER)
   
where ORDER = p+1 (order of error estimate)
```

**Safety factor**: Typically 0.9 to prevent oscillation

**Growth/shrink limits**: Prevent h from changing too rapidly

### Symplectic Integrators

**For Hamiltonian systems**: Preserve phase space volume, energy (approximately)

**Leapfrog** (also called Störmer-Verlet):
```
Split: y = [q, p]  (position, momentum)

v_{n+1/2} = v_n + (h/2)·a(q_n)     # Half kick
q_{n+1} = q_n + h·v_{n+1/2}         # Full drift
v_{n+1} = v_{n+1/2} + (h/2)·a(q_{n+1})  # Half kick
```

**Properties**:
- Time-reversible
- Symplectic (preserves phase space structure)
- Energy error bounded (doesn't grow)
- 2nd-order accuracy

**Velocity Verlet** (equivalent formulation):
```
q_{n+1} = q_n + h·v_n + (h²/2)·a(q_n)
v_{n+1} = v_n + (h/2)·(a(q_n) + a(q_{n+1}))
```

---

## Classes and Functions

### IODESystem - Interface for ODE Systems

**Purpose**: Define dy/dt = f(t,y)

```cpp
class IODESystem {
public:
    virtual int getDim() const = 0;  // System dimension
    
    // Compute derivatives: dydt = f(t, y)
    virtual void derivs(Real t, const Vector<Real>& y, Vector<Real>& dydt) const = 0;
};
```

**Implementation Example**:
```cpp
class LorenzSystem : public IODESystem {
    Real sigma, rho, beta;
public:
    LorenzSystem(Real s, Real r, Real b) : sigma(s), rho(r), beta(b) {}
    
    int getDim() const override { return 3; }
    
    void derivs(Real t, const Vector<Real>& y, Vector<Real>& dydt) const override {
        dydt[0] = sigma * (y[1] - y[0]);
        dydt[1] = y[0] * (rho - y[2]) - y[1];
        dydt[2] = y[0] * y[1] - beta * y[2];
    }
};
```

---

### IODESystemStepCalculator - Interface for Step Calculators

**Purpose**: Single step computation

```cpp
class IODESystemStepCalculator {
public:
    virtual void calcStep(
        const IODESystem& sys,
        const Real t,                    // Current time
        const Vector<Real>& x_start,     // Current state
        const Vector<Real>& dxdt,        // Current derivatives
        const Real h,                    // Step size
        Vector<Real>& x_out,             // OUTPUT: New state
        Vector<Real>& x_err_out          // OUTPUT: Error estimate (if available)
    ) const = 0;
};
```

---

### Step Calculators (Fixed-Step Methods)

#### EulerStep_Calculator

**Order**: 1  
**Evaluations**: 1 per step

```cpp
class EulerStep_Calculator : public IODESystemStepCalculator {
    void calcStep(...) const override {
        for (int i = 0; i < n; i++)
            x_out[i] = x_start[i] + h * dxdt[i];
    }
};
```

**Formula**: y_{n+1} = y_n + h·f(t_n, y_n)

**When to use**: Teaching, debugging, very simple problems  
**Accuracy**: Poor (O(h))  
**Stability**: Bad

#### Midpoint_StepCalculator

**Order**: 2  
**Evaluations**: 2 per step

```cpp
class Midpoint_StepCalculator : public IODESystemStepCalculator {
    void calcStep(...) const override {
        Vector<Real> x_mid(n);
        
        // Step to midpoint
        for (int i = 0; i < n; i++)
            x_mid[i] = x_start[i] + 0.5 * h * dxdt[i];
        
        // Evaluate at midpoint, use for full step
        sys.derivs(t + 0.5*h, x_mid, x_out);
    }
};
```

**Formula**: 
```
y_mid = y_n + (h/2)·f(t_n, y_n)
y_{n+1} = y_n + h·f(t_n + h/2, y_mid)
```

**When to use**: Simple problems, better than Euler  
**Accuracy**: O(h²)

#### RungeKutta4_StepCalculator (RK4)

**Order**: 4 (classic workhorse!)  
**Evaluations**: 4 per step

```cpp
class RungeKutta4_StepCalculator : public IODESystemStepCalculator {
    void calcStep(...) const override {
        Vector<Real> k1, k2, k3, k4, x_temp;
        
        k1 = dxdt;
        
        // k2 at t + h/2
        x_temp = x_start + (h/2)*k1;
        sys.derivs(t + h/2, x_temp, k2);
        
        // k3 at t + h/2
        x_temp = x_start + (h/2)*k2;
        sys.derivs(t + h/2, x_temp, k3);
        
        // k4 at t + h
        x_temp = x_start + h*k3;
        sys.derivs(t + h, x_temp, k4);
        
        // Weighted combination
        x_out = x_start + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    }
};
```

**Formula**:
```
k₁ = f(t, y)
k₂ = f(t + h/2, y + h·k₁/2)
k₃ = f(t + h/2, y + h·k₂/2)
k₄ = f(t + h, y + h·k₃)
y_{n+1} = y_n + (h/6)(k₁ + 2k₂ + 2k₃ + k₄)
```

**When to use**: General-purpose fixed-step, standard method  
**Accuracy**: O(h⁴) - excellent for fixed step  
**Stability**: Good for non-stiff

#### VelocityVerlet_StepCalculator

**Order**: 2  
**Evaluations**: 2 per step  
**For**: Hamiltonian systems (position/velocity split)

```cpp
class VelocityVerlet_StepCalculator : public IODESystemStepCalculator {
    void calcStep(...) const override {
        int half_n = n / 2;
        
        // State: [positions(0..half_n-1), velocities(half_n..n-1)]
        // Derivs: [velocities(0..half_n-1), accelerations(half_n..n-1)]
        
        // 1. Update positions
        for (int i = 0; i < half_n; ++i)
            x_out[i] = x_start[i] + x_start[half_n+i]*h + 0.5*dxdt[half_n+i]*h*h;
        
        // 2. Compute new accelerations
        Vector<Real> dxdt_new(n);
        sys.derivs(t + h, x_out, dxdt_new);
        
        // 3. Update velocities (average of old and new accelerations)
        for (int i = 0; i < half_n; ++i)
            x_out[half_n+i] = x_start[half_n+i] + 0.5*(dxdt[half_n+i] + dxdt_new[half_n+i])*h;
    }
};
```

**Formula**:
```
q_{n+1} = q_n + h·v_n + (h²/2)·a_n
v_{n+1} = v_n + (h/2)·(a_n + a_{n+1})
```

**Properties**:
- Symplectic
- Energy-conserving (bounded error)
- Time-reversible

**When to use**: Molecular dynamics, planetary motion, Hamiltonian mechanics

#### Leapfrog_StepCalculator

**Order**: 2  
**Evaluations**: 2 per step  
**For**: Hamiltonian systems

```cpp
class Leapfrog_StepCalculator : public IODESystemStepCalculator {
    void calcStep(...) const override {
        // 1. Half-kick: v_{n+1/2} = v_n + (h/2)*a_n
        // 2. Full-drift: q_{n+1} = q_n + h*v_{n+1/2}
        // 3. Compute a_{n+1}
        // 4. Half-kick: v_{n+1} = v_{n+1/2} + (h/2)*a_{n+1}
    }
};
```

**When to use**: Same as Velocity Verlet (equivalent for constant timestep)

---

### Adaptive Step Calculators (Embedded Methods)

#### RK5_CashKarp_Calculator

**Order**: 5(4) embedded pair  
**Evaluations**: 6 per step

**Butcher Tableau**: Cash-Karp coefficients

```cpp
class RK5_CashKarp_Calculator : public IODESystemStepCalculator {
    void calcStep(...) const override {
        // Compute 6 intermediate k values
        // 5th-order solution: weighted combination
        // 4th-order solution: different weights
        // Error = ||5th - 4th||
    }
};
```

**When to use**: Good general-purpose adaptive method  
**Properties**:
- Robust error control
- Efficient
- Not FSAL (need 6 evals even if step rejected)

#### DormandPrince5_StepCalculator

**Order**: 5(4) embedded pair  
**Evaluations**: 7 per step (but FSAL → effectively 6)

```cpp
class DormandPrince5_StepCalculator : public IODESystemStepCalculator {
    void calcStep(...) const override {
        // Compute 7 intermediate k values
        // FSAL property: k₇ becomes k₁ of next step
        // 5th-order solution for y_{n+1}
        // 4th-order for error estimate
    }
};
```

**FSAL** (First Same As Last):
- Last evaluation k₇ of current step = k₁ of next step
- Saves 1 function evaluation per step
- **Most efficient general-purpose method**

**When to use**: **DEFAULT CHOICE** for non-stiff problems  
**Industry standard**: MATLAB ode45, SciPy solve_ivp default

#### DormandPrince8_StepCalculator

**Order**: 8(7) embedded pair  
**Evaluations**: 13 per step

```cpp
class DormandPrince8_StepCalculator : public IODESystemStepCalculator {
    void calcStep(...) const override {
        // 13 intermediate k values
        // 8th-order solution
        // 7th-order error estimate
    }
};
```

**When to use**: High-precision requirements (ε < 10⁻⁸)  
**Trade-off**: More evaluations per step, but can use larger h

---

### Steppers (Adaptive Control Wrappers)

#### StepperBase

**Abstract base** for adaptive steppers

```cpp
class StepperBase {
protected:
    const IODESystem& _sys;
    Real& _t;              // Reference to solver's time
    Vector<Real>& _x;      // Reference to solver's state
    Vector<Real>& _dxdt;   // Reference to solver's derivatives
    
    Real _hDone;           // Actual step taken
    Real _hNext;           // Recommended next step
    
public:
    virtual void doStep(Real htry, Real eps) = 0;
};
```

**Responsibilities**:
- Call step calculator
- Check error
- Accept/reject step
- Adjust timestep

#### RK5_CashKarp_Stepper

```cpp
class RK5_CashKarp_Stepper : public StepperBase {
    RK5_CashKarp_Calculator _stepCalc;
    
public:
    void doStep(Real htry, Real eps) override {
        while (true) {
            _stepCalc.calcStep(_sys, _t, _x, _dxdt, htry, _xout, _xerr);
            
            Real errmax = max_i(|xerr[i]| / scale[i]) / eps;
            
            if (errmax <= 1.0) {
                // ACCEPT
                _hDone = htry;
                _hNext = SAFETY * htry * pow(errmax, PGROW);
                _x = _xout;
                _t += htry;
                return;
            }
            
            // REJECT - retry with smaller step
            htry = SAFETY * htry * pow(errmax, PSHRNK);
        }
    }
};
```

**Constants**:
- SAFETY = 0.9 (prevent oscillation)
- PGROW = -0.2 (exponent for growth)
- PSHRNK = -0.25 (exponent for shrinkage)

---

### Solvers

#### ODESystemFixedStepSolver

**For**: Fixed timestep integration

```cpp
class ODESystemFixedStepSolver {
    const IODESystem& _sys;
    const IODESystemStepCalculator& _stepCalc;
    
public:
    ODESystemSolution integrate(
        const Vector<Real>& initCond,
        Real t1, Real t2,
        int numSteps);
};
```

**Algorithm**:
```
h = (t2 - t1) / numSteps
for k = 1 to numSteps:
    stepCalc.calcStep(sys, t, x, dxdt, h, x_out, x_err)
    t += h
    x = x_out
    store (t, x) in solution
```

**Returns**: `ODESystemSolution` with exactly numSteps+1 points

#### ODESystemSolver<Stepper>

**For**: Adaptive timestep integration

```cpp
template<class Stepper>
class ODESystemSolver {
    const IODESystem& _sys;
    Stepper _stepper;
    
    Real _curr_t;
    Vector<Real> _curr_x;
    Vector<Real> _curr_dxdt;
    
public:
    ODESystemSolution integrate(
        const Vector<Real>& initCond,
        Real t1, Real t2,
        Real minSaveInterval,  // Minimum spacing between saved points
        Real eps,              // Error tolerance
        Real h1,               // Initial step size
        Real hmin = 0);        // Minimum allowed step
};
```

**Algorithm**:
```
t = t1, x = initCond
while (t < t2):
    stepper.doStep(htry, eps)
    if (t - lastSave >= minSaveInterval):
        store (t, x) in solution
    htry = stepper.hNext()
```

**Features**:
- Automatic error control
- Dense output (interpolated points at regular intervals)
- Step rejection for accuracy
- Efficient (large steps when possible)

**Typical usage**:
```cpp
LorenzSystem sys(10.0, 28.0, 8.0/3.0);
ODESystemSolver<RK5_CashKarp_Stepper> solver(sys);

Vector<Real> y0 = {1.0, 1.0, 1.0};
Real t1 = 0.0, t2 = 50.0;
Real minSave = 0.01;  // Save every 0.01 units
Real eps = 1e-6;      // 10^-6 accuracy
Real h1 = 0.01;       // Initial step

auto sol = solver.integrate(y0, t1, t2, minSave, eps, h1);
```

#### ODESystemLeapfrogSolver

**Specialized** for Hamiltonian systems

```cpp
class ODESystemLeapfrogSolver {
    IODESystem& _sys;
    
public:
    ODESystemSolution integrate(
        const Vector<Real>& initCond,  // [positions, velocities]
        Real t1, Real t2,
        int numSteps);
};
```

**Fixed-step**, symplectic Leapfrog integration

---

## Solution Output

### ODESystemSolution

```cpp
class ODESystemSolution {
public:
    Vector<Real> _xval;     // Time points (length: nSave)
    MatrixNM<Real> _yval;   // Solution (nSave × dim)
    
    Real t1() const;        // Start time
    Real t2() const;        // End time
    int dim() const;        // System dimension
    int numSavedSteps() const;  // Number of saved points
    
    void fillValues(int step, Real t, const Vector<Real>& y);
};
```

**Access**:
```cpp
Real t_i = sol._xval[i];
Real y_j_at_i = sol._yval[i][j];  // Component j at time i
```

---

## Examples

### Example 1: Simple Harmonic Oscillator (Fixed RK4)

```cpp
#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"

// d²x/dt² = -ω²x  →  dy/dt = [v, -ω²y]
class SHO : public IODESystem {
    Real omega;
public:
    SHO(Real w) : omega(w) {}
    int getDim() const override { return 2; }
    
    void derivs(Real t, const Vector<Real>& y, Vector<Real>& dydt) const override {
        dydt[0] = y[1];              // dx/dt = v
        dydt[1] = -omega*omega*y[0]; // dv/dt = -ω²x
    }
};

void Example1() {
    SHO sys(1.0);  // ω = 1
    RungeKutta4_StepCalculator rk4;
    ODESystemFixedStepSolver solver(sys, rk4);
    
    Vector<Real> y0 = {1.0, 0.0};  // x=1, v=0
    auto sol = solver.integrate(y0, 0.0, 10.0, 1000);
    
    // Exact: x(t) = cos(t), v(t) = -sin(t)
    for (int i = 0; i < sol.numSavedSteps(); ++i) {
        Real t = sol._xval[i];
        Real x = sol._yval[i][0];
        Real x_exact = std::cos(t);
        Real error = std::abs(x - x_exact);
        
        if (i % 100 == 0)
            std::cout << "t=" << t << ", x=" << x << ", error=" << error << "\n";
    }
}
```

### Example 2: Lorenz System (Adaptive Dormand-Prince)

Chaotic system - needs adaptive stepping!

```cpp
class LorenzSystem : public IODESystem {
    Real sigma, rho, beta;
public:
    LorenzSystem(Real s, Real r, Real b) : sigma(s), rho(r), beta(b) {}
    int getDim() const override { return 3; }
    
    void derivs(Real t, const Vector<Real>& y, Vector<Real>& dydt) const override {
        dydt[0] = sigma * (y[1] - y[0]);
        dydt[1] = y[0] * (rho - y[2]) - y[1];
        dydt[2] = y[0] * y[1] - beta * y[2];
    }
};

void Example2() {
    LorenzSystem sys(10.0, 28.0, 8.0/3.0);  // Classic parameters
    ODESystemSolver<RK5_CashKarp_Stepper> solver(sys);
    
    Vector<Real> y0 = {1.0, 1.0, 1.0};
    Real t1 = 0.0, t2 = 50.0;
    Real minSave = 0.01;  // Dense output every 0.01
    Real eps = 1e-6;      // High accuracy
    Real h1 = 0.01;
    
    auto sol = solver.integrate(y0, t1, t2, minSave, eps, h1);
    
    std::cout << "Solved in " << sol.numSavedSteps() << " saved steps\n";
    
    // Visualize (write to file for plotting)
    std::ofstream out("lorenz.txt");
    for (int i = 0; i < sol.numSavedSteps(); ++i) {
        out << sol._xval[i] << " "
            << sol._yval[i][0] << " "
            << sol._yval[i][1] << " "
            << sol._yval[i][2] << "\n";
    }
}
```

### Example 3: Comparison of Methods

Same problem, different solvers:

```cpp
void Example3() {
    SHO sys(2.0);  // ω = 2
    Vector<Real> y0 = {1.0, 0.0};
    Real t1 = 0.0, t2 = 5.0;
    
    // Method 1: Euler (terrible)
    {
        EulerStep_Calculator euler;
        ODESystemFixedStepSolver solver(sys, euler);
        auto sol = solver.integrate(y0, t1, t2, 1000);
        std::cout << "Euler final x: " << sol._yval[1000][0] << " (exact: " << std::cos(2*5) << ")\n";
    }
    
    // Method 2: RK4 (good)
    {
        RungeKutta4_StepCalculator rk4;
        ODESystemFixedStepSolver solver(sys, rk4);
        auto sol = solver.integrate(y0, t1, t2, 100);  // 10× fewer steps!
        std::cout << "RK4 final x: " << sol._yval[100][0] << " (exact: " << std::cos(2*5) << ")\n";
    }
    
    // Method 3: Adaptive (best)
    {
        ODESystemSolver<RK5_CashKarp_Stepper> solver(sys);
        auto sol = solver.integrate(y0, t1, t2, 0.05, 1e-8, 0.1);
        int last = sol.numSavedSteps() - 1;
        std::cout << "Adaptive final x: " << sol._yval[last][0] << " (exact: " << std::cos(2*5) << ")\n";
        std::cout << "Used " << sol.numSavedSteps() << " saved points\n";
    }
}
```

### Example 4: Planetary Motion (Leapfrog)

Two-body problem - energy conservation critical!

```cpp
class TwoBodyProblem : public IODESystem {
    Real G, M;  // Gravitational constant, central mass
public:
    TwoBodyProblem(Real g, Real m) : G(g), M(m) {}
    int getDim() const override { return 4; }  // [x, y, vx, vy]
    
    void derivs(Real t, const Vector<Real>& y, Vector<Real>& dydt) const override {
        Real x = y[0], y_pos = y[1];
        Real vx = y[2], vy = y[3];
        
        Real r = std::sqrt(x*x + y_pos*y_pos);
        Real r3 = r*r*r;
        
        dydt[0] = vx;
        dydt[1] = vy;
        dydt[2] = -G*M*x / r3;  // ax
        dydt[3] = -G*M*y_pos / r3;  // ay
    }
};

void Example4() {
    TwoBodyProblem sys(1.0, 1.0);  // G=1, M=1
    
    // Circular orbit: v = sqrt(GM/r)
    Real r0 = 1.0;
    Real v0 = std::sqrt(1.0 / r0);
    Vector<Real> y0 = {r0, 0.0, 0.0, v0};
    
    // Use Leapfrog for energy conservation
    VelocityVerlet_StepCalculator verlet;
    ODESystemFixedStepSolver solver(sys, verlet);
    
    Real T = 2 * Constants::PI;  // One orbital period
    auto sol = solver.integrate(y0, 0.0, 10*T, 10000);
    
    // Check energy conservation
    for (int i = 0; i < sol.numSavedSteps(); i += 1000) {
        Real x = sol._yval[i][0], y = sol._yval[i][1];
        Real vx = sol._yval[i][2], vy = sol._yval[i][3];
        
        Real KE = 0.5 * (vx*vx + vy*vy);
        Real PE = -1.0 / std::sqrt(x*x + y*y);
        Real E = KE + PE;
        
        std::cout << "t=" << sol._xval[i] << ", E=" << E << "\n";
    }
    // Energy should be constant (E = -0.5 for circular orbit at r=1)
}
```

### Example 5: Stiff Problem (Chemical Kinetics)

Robertson's problem - classic stiff ODE

```cpp
class RobertsonProblem : public IODESystem {
public:
    int getDim() const override { return 3; }
    
    void derivs(Real t, const Vector<Real>& y, Vector<Real>& dydt) const override {
        // dy1/dt = -0.04*y1 + 1e4*y2*y3
        // dy2/dt = 0.04*y1 - 1e4*y2*y3 - 3e7*y2^2
        // dy3/dt = 3e7*y2^2
        
        dydt[0] = -0.04*y[0] + 1e4*y[1]*y[2];
        dydt[1] = 0.04*y[0] - 1e4*y[1]*y[2] - 3e7*y[1]*y[1];
        dydt[2] = 3e7*y[1]*y[1];
    }
};

void Example5() {
    RobertsonProblem sys;
    Vector<Real> y0 = {1.0, 0.0, 0.0};
    
    // Stiff! Need small steps even with adaptive method
    ODESystemSolver<RK5_CashKarp_Stepper> solver(sys);
    
    Real eps = 1e-4;  // Looser tolerance for stiff problem
    Real h1 = 1e-6;   // Very small initial step
    Real hmin = 1e-10; // Allow tiny steps
    
    auto sol = solver.integrate(y0, 0.0, 1e5, 1.0, eps, h1, hmin);
    
    std::cout << "Solved stiff problem in " << sol.numSavedSteps() << " steps\n";
    
    // Check conservation: y1 + y2 + y3 = 1
    int last = sol.numSavedSteps() - 1;
    Real sum = sol._yval[last][0] + sol._yval[last][1] + sol._yval[last][2];
    std::cout << "Conservation check: sum = " << sum << " (should be 1.0)\n";
}
```

### Example 6: High Precision (Dormand-Prince 8)

When you need 10+ digits:

```cpp
void Example6() {
    SHO sys(1.0);
    Vector<Real> y0 = {1.0, 0.0};
    
    ODESystemSolver<DormandPrince8_Stepper> solver(sys);
    
    Real eps = 1e-12;  // Very tight tolerance
    auto sol = solver.integrate(y0, 0.0, 10.0, 0.1, eps, 0.01);
    
    // Check final accuracy
    int last = sol.numSavedSteps() - 1;
    Real t_final = sol._xval[last];
    Real x_computed = sol._yval[last][0];
    Real x_exact = std::cos(t_final);
    
    std::cout << std::setprecision(15);
    std::cout << "Computed: " << x_computed << "\n";
    std::cout << "Exact:    " << x_exact << "\n";
    std::cout << "Error:    " << std::abs(x_computed - x_exact) << "\n";
    // Typical: error ~ 1e-12
}
```

---

## Integration with Other Modules

### With TestBeds

**Predefined systems** for testing:

```cpp
namespace TestBeds {
    class LorenzSystemODE : public IODESystem { /* ... */ };
    class DoublePendulum : public IODESystem { /* ... */ };
    class VanDerPolOscillator : public IODESystem { /* ... */ };
}

// Usage:
TestBeds::LorenzSystemODE sys(10.0, 28.0, 8.0/3.0);
ODESystemSolver<RK5_CashKarp_Stepper> solver(sys);
// ...
```

### With Visualizers

**Export solutions** for plotting:

```cpp
// Visualizer namespace provides helpers
Visualizer::VisualizeODESysSolAsMultiFunc(
    sol, 
    "Lorenz System Components vs Time", 
    "lorenz_components.txt");

Visualizer::VisualizeODESysSolAsParamCurve3(
    sol,
    "Lorenz Attractor 3D",
    "lorenz_attractor_3d.txt");
```

**Output format**: Space-separated text files for gnuplot, Python, etc.

---

## Best Practices

### Choosing a Method

**Decision Tree**:

```
Is your problem Hamiltonian (energy-conserving)?
├─ YES → Use Leapfrog or VelocityVerlet
└─ NO → Continue

Do you know the problem is non-stiff?
├─ YES → Use Dormand-Prince 5 (adaptive)
└─ UNSURE → Try DP5, monitor step sizes

Do you need very high precision (< 10⁻⁸)?
├─ YES → Use Dormand-Prince 8
└─ NO → Dormand-Prince 5 is fine

Is computational cost critical?
├─ YES → Use RK4 fixed-step (if accuracy acceptable)
└─ NO → Use adaptive (safer)

Is problem known to be stiff?
└─ Use small timesteps, looser tolerance
   (Better: use implicit methods when available)
```

### Error Tolerance Selection

**Rule of Thumb**:
```
For relative tolerance ε:
  Expect accuracy ~ ε × ||y||

ε = 1e-3:  3 significant digits
ε = 1e-6:  6 significant digits
ε = 1e-9:  9 significant digits
```

**Trade-offs**:
- Smaller ε → More steps, higher cost
- Larger ε → Faster, less accurate

### Initial Step Size

**Good heuristics**:
```cpp
Real h1 = 0.01 * (t2 - t1);  // 1% of integration interval

// Or estimate from derivative:
Vector<Real> dydt(n);
sys.derivs(t1, y0, dydt);
Real scale = maxnorm(y0);  // or 1.0 if y0 small
Real deriv_scale = maxnorm(dydt);
Real h1 = eps * scale / deriv_scale;
```

**Note**: Adaptive methods adjust automatically, so initial h1 not critical

### Monitoring Warnings

**Watch for**:
```cpp
// If many step rejections:
if (num_steps_taken >> expected_steps)
    std::cerr << "Warning: Many step rejections - problem may be stiff!\n";

// If hitting minimum step size:
if (h < hmin)
    std::cerr << "Warning: Hit minimum step size - may lose accuracy!\n";

// Non-finite values:
if (std::isnan(y) || std::isinf(y))
    throw ODESolverError("Non-finite value detected");
```

### Common Pitfalls

❌ **Pitfall 1**: Using Euler for production
```cpp
// DON'T:
EulerStep_Calculator euler;
solver(sys, euler).integrate(...);  // Terrible accuracy!

// DO:
RungeKutta4_StepCalculator rk4;  // At minimum
// Or better: adaptive DP5
```

❌ **Pitfall 2**: Wrong method for problem type
```cpp
// DON'T (Hamiltonian system):
ODESystemSolver<RK5_CashKarp_Stepper> solver(planetary_system);
// Energy will drift!

// DO:
VelocityVerlet_StepCalculator verlet;
ODESystemFixedStepSolver solver(planetary_system, verlet);
// Energy conserved (bounded error)
```

❌ **Pitfall 3**: Too tight tolerance
```cpp
// DON'T:
Real eps = 1e-15;  // Beyond machine precision!
solver.integrate(..., eps, ...);

// DO:
Real eps = 1e-10;  // Reasonable for double precision
```

❌ **Pitfall 4**: Not converting second-order equations
```cpp
// DON'T: Try to solve d²y/dt² = f(y) directly

// DO: Convert to first-order system:
// y₁ = y, y₂ = dy/dt
// dy₁/dt = y₂
// dy₂/dt = f(y₁)
```

✅ **Best Practice Pattern**:
```cpp
// 1. Define system clearly
class MySystem : public IODESystem {
    // ...
};

// 2. Choose appropriate method
MySystem sys(...);
ODESystemSolver<DormandPrince5_Stepper> solver(sys);  // Default choice

// 3. Set reasonable parameters
Real eps = 1e-6;      // 6-digit accuracy
Real h1 = 0.01;       // Initial step
Real minSave = 0.01;  // Dense output

// 4. Integrate
auto sol = solver.integrate(y0, t1, t2, minSave, eps, h1);

// 5. Validate
if (sol.numSavedSteps() == 0)
    throw std::runtime_error("Integration failed!");

// 6. Check conservation laws if applicable
// (energy, momentum, etc.)
```

---

## Performance Considerations

### Computational Cost

| Method | Evals/Step | Order | Cost for ε accuracy |
|--------|-----------|-------|---------------------|
| Euler | 1 | 1 | O((t₂-t₁)/ε) |
| RK4 | 4 | 4 | O((t₂-t₁)/ε^(1/4)) |
| DP5 | 6* | 5 | O((t₂-t₁)/ε^(1/5)) |
| DP8 | 13 | 8 | O((t₂-t₁)/ε^(1/8)) |

*FSAL: effectively 6 after first step

**Example**: For ε = 10⁻⁶, interval [0,10]
- Euler: ~10⁷ function evaluations
- RK4: ~40,000 evaluations
- DP5 adaptive: ~6,000 evaluations
- DP8 adaptive: ~2,000 evaluations (but 13 evals/step)

**Conclusion**: DP5 often optimal balance

### When to Use Fixed vs Adaptive

**Fixed-step (RK4)**:
- Simple, predictable problems
- Hamiltonian systems (with Verlet/Leapfrog)
- Real-time applications (predictable timing)
- Benchmarking
- Teaching

**Adaptive (DP5/DP8)**:
- Unknown behavior
- Accuracy critical
- Efficiency matters
- Stiff regions possible
- Production code

### Function Evaluation Optimization

**Most cost** is in `sys.derivs()` calls - optimize there!

```cpp
// Slow:
void derivs(Real t, const Vector<Real>& y, Vector<Real>& dydt) const {
    Real expensive = compute_expensive_thing(y);  // Called for EACH component!
    dydt[0] = expensive * y[1];
    dydt[1] = expensive * y[0];
}

// Fast:
void derivs(Real t, const Vector<Real>& y, Vector<Real>& dydt) const {
    Real expensive = compute_expensive_thing(y);  // Called ONCE
    dydt[0] = expensive * y[1];
    dydt[1] = expensive * y[0];
}
```

---

## Advanced Topics

### Dense Output

**Problem**: Adaptive methods produce **irregular** time points

**Solution**: Interpolation for dense output at regular intervals

**Implementation**: `minSaveInterval` parameter
```cpp
solver.integrate(y0, 0, 10, 
                 0.01,  // Save every 0.01 (dense output)
                 1e-6, 0.01);
```

**How it works**: 
- Solver takes adaptive steps (irregular)
- Whenever t - t_lastSave > minSaveInterval, interpolate and save
- Result: Regular output grid for visualization/post-processing

### Discontinuous Right-Hand Side

**Problem**: Discontinuities in f(t,y) (collisions, switches)

**Example**: Bouncing ball (velocity reverses at ground)

**Solution**: 
1. Detect event (y = 0)
2. Stop integration at event
3. Modify state (velocity reversal)
4. Restart integration

```cpp
// Detect ground collision
if (y[0] <= 0 && v[1] < 0) {
    y[1] = -0.9 * y[1];  // Coefficient of restitution 0.9
}
```

### Jacobian for Stiff Problems

**Stiff solvers** (implicit methods, not yet implemented) need Jacobian:

**J_ij = ∂f_i/∂y_j**

**Can provide** as additional interface:
```cpp
class StiffSystem : public IODESystem {
    void jacobian(Real t, const Vector<Real>& y, MatrixNM<Real>& J) const {
        // Compute ∂f/∂y
    }
};
```

**Future**: Rosenbrock, BDF methods for stiff problems

---

## Summary

### Method Quick Reference

**For Hamiltonian Systems**:
- Leapfrog or Velocity Verlet (symplectic, energy-conserving)

**For General Non-Stiff**:
- **First choice**: Dormand-Prince 5 (adaptive)
- High precision: Dormand-Prince 8
- Simple/teaching: RK4 fixed-step

**For Stiff Systems**:
- Small steps with DP5
- Better: Implicit methods (future work)

### Key Takeaways

1. ✅ **Separate concerns**: Calculator (math) vs Stepper (control) vs Solver (orchestration)
2. ✅ **Adaptive >> Fixed** for most production code
3. ✅ **DP5 is industry standard** for non-stiff ODEs
4. ✅ **Symplectic methods preserve structure** (Hamiltonian systems)
5. ✅ **Error tolerance**: ε ~ 10⁻⁶ good default
6. ✅ **Watch for stiffness**: Many rejected steps = problem is stiff
7. ✅ **Validate**: Check conservation laws, compare with known solutions

### References

- **Numerical Recipes** (Press et al.): Chapter 16 - Integration of ODEs
- **Hairer, Nørsett, Wanner**: *Solving Ordinary Differential Equations I* (Non-Stiff)
- **Hairer, Wanner**: *Solving Ordinary Differential Equations II* (Stiff)
- **Dormand & Prince (1980)**: Original DP5(4) paper
- **Cash & Karp (1990)**: RK5(4) embedded method

---

## Runnable Examples

Working code examples are available in `src/docs_demos/docs_demo_ode_solvers.cpp`:

### IODESystem Interface
- `Docs_Demo_IODESystem_Interface()` - Creating custom ODE systems

### ODESystem Creation Patterns
- `Docs_Demo_ODESystem_Creation()` - Different ways to define ODE systems

### Step Calculators
- `Docs_Demo_Step_Calculators()` - Euler, Midpoint, RK4, Leapfrog comparisons

### Adaptive Integration
- `Docs_Demo_Adaptive_Integrators()` - Cash-Karp, Dormand-Prince 5/8

### Legacy Interface
- `Docs_Demo_Legacy_Solver_Interface()` - Backward-compatible ODESystemSolver

### Solution Post-Processing
- `Docs_Demo_ODESystemSolution_Basics()` - Accessing solution data
- `Docs_Demo_ODESystemSolution_Interpolation()` - Interpolating between saved points
- `Docs_Demo_ODESystemSolution_ParametricCurve()` - Converting to parametric curves

### Physics Examples
- `Docs_Demo_Physics_SHO()` - Simple harmonic oscillator
- `Docs_Demo_Physics_Damped_Oscillator()` - Damped harmonic oscillator
- `Docs_Demo_Physics_Lorenz()` - Lorenz attractor (chaotic system)
- `Docs_Demo_Physics_Predator_Prey()` - Lotka-Volterra equations

### Classic Demo
- `Docs_Demo_ODE_solvers_RungeKutta_4th_order()` - RK4 and adaptive solver comparison

---

**Part of MinimalMathLibrary** - `mml/algorithms/ODESystem*.h`
