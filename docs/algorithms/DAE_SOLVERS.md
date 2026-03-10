# DAE Solvers

Differential-Algebraic Equation (DAE) solvers for semi-explicit index-1 systems.

## Overview

DAE systems combine differential equations with algebraic constraints:

```
dx/dt = f(t, x, y)   (differential equations)
0 = g(t, x, y)       (algebraic constraints)
```

Where:
- **x** = differential variables (have time derivatives)
- **y** = algebraic variables (determined instantaneously by constraints)
- **f** = differential right-hand side
- **g** = algebraic constraint equations

Common examples include:
- Constrained mechanical systems (pendulum, robotics)
- Electrical circuits with ideal components
- Chemical reaction networks with conservation laws
- Incompressible fluid dynamics

## Index Concept

The **index** of a DAE measures the "distance" from an ODE:
- **Index-0**: Pure ODE (no constraints)
- **Index-1**: Constraints involve only algebraic variables directly
- **Index-2+**: Higher-order DAEs requiring more derivatives to reduce

The MML solvers handle **index-1 systems** where ∂g/∂y is nonsingular.

## Solver Methods

### Backward Euler (1st Order)

```cpp
#include "mml/algorithms/DAESolvers.h"

// Define your DAE system
class MyDAE : public IODESystemDAEWithJacobian {
    int getDiffDim() const override { return 2; }
    int getAlgDim() const override { return 1; }
    
    void diffEqs(Real t, const Vector<Real>& x, const Vector<Real>& y,
                 Vector<Real>& dxdt) const override {
        // ... your differential equations
    }
    
    void algConstraints(Real t, const Vector<Real>& x, const Vector<Real>& y,
                        Vector<Real>& g) const override {
        // ... your constraint equations
    }
    
    // ... Jacobian methods
};

// Solve
MyDAE system;
Vector<Real> x0(2), y0(1);  // Initial conditions
x0[0] = 1.0; x0[1] = 0.0;
y0[0] = 0.0;

DAESolverConfig config;
config.step_size = 0.01;

DAESolverResult result = SolveDAEBackwardEuler(system, 0.0, x0, y0, 10.0, config);
```

Properties:
- **Order**: 1 (error ∝ h)
- **Stability**: A-stable, L-stable
- **Best for**: Stiff problems, simple baseline

### BDF2 (2nd Order)

```cpp
DAESolverResult result = SolveDAEBDF2(system, 0.0, x0, y0, 10.0, config);
```

Properties:
- **Order**: 2 (error ∝ h²)
- **Stability**: A-stable
- **Best for**: Higher accuracy with moderate stiffness

The BDF2 method uses Backward Euler for the first step to bootstrap.

## Configuration

```cpp
struct DAESolverConfig {
    Real step_size = 0.01;          // Fixed step size
    int max_newton_iter = 20;       // Newton iterations per step
    Real newton_tol = 1e-8;         // Newton convergence tolerance
    Real constraint_tol = 1e-10;    // Constraint satisfaction tolerance
    int max_steps = 100000;         // Maximum integration steps
};

// Preset configurations
DAESolverConfig highPrec = DAESolverConfig::HighPrecision();
DAESolverConfig fast = DAESolverConfig::Fast();
```

## Result Structure

```cpp
struct DAESolverResult {
    DAESolution solution;           // Solution trajectory
    
    std::string algorithm_name;     // "DAEBackwardEuler" or "DAEBDF2"
    AlgorithmStatus status;         // Success, NumericalInstability, etc.
    std::string error_message;
    double elapsed_time_ms;
    
    int total_steps;
    int newton_iterations;
    int jacobian_evaluations;
    Real max_constraint_violation;
};
```

## Consistent Initial Conditions

DAE solvers require **consistent initial conditions** where g(t₀, x₀, y₀) = 0.

### Verify Consistency

```cpp
bool ok = VerifyConsistentIC(system, t0, x0, y0, 1e-10);
```

### Compute Consistent y₀

Given x₀, solve for y₀ such that constraints are satisfied:

```cpp
Vector<Real> y0(1);
y0[0] = 0.5;  // Initial guess

bool success = ComputeConsistentIC(system, t0, x0, y0, /*max_iter=*/20, /*tol=*/1e-10);
// y0 now satisfies g(t0, x0, y0) = 0
```

## DAE System Interface

### Basic Interface

```cpp
class IODESystemDAE {
public:
    virtual int getDiffDim() const = 0;  // Size of x
    virtual int getAlgDim() const = 0;   // Size of y
    
    // dx/dt = f(t, x, y)
    virtual void diffEqs(Real t, const Vector<Real>& x, const Vector<Real>& y,
                         Vector<Real>& dxdt) const = 0;
    
    // 0 = g(t, x, y)
    virtual void algConstraints(Real t, const Vector<Real>& x, const Vector<Real>& y,
                                 Vector<Real>& g) const = 0;
};
```

### With Jacobians (Required for Implicit Solvers)

```cpp
class IODESystemDAEWithJacobian : public IODESystemDAE {
public:
    virtual void jacobian_fx(Real t, const Vector<Real>& x, const Vector<Real>& y,
                             Matrix<Real>& df_dx) const = 0;
    
    virtual void jacobian_fy(Real t, const Vector<Real>& x, const Vector<Real>& y,
                             Matrix<Real>& df_dy) const = 0;
    
    virtual void jacobian_gx(Real t, const Vector<Real>& x, const Vector<Real>& y,
                             Matrix<Real>& dg_dx) const = 0;
    
    virtual void jacobian_gy(Real t, const Vector<Real>& x, const Vector<Real>& y,
                             Matrix<Real>& dg_dy) const = 0;
    
    // Convenience method for all four
    virtual void allJacobians(Real t, const Vector<Real>& x, const Vector<Real>& y,
                              Matrix<Real>& df_dx, Matrix<Real>& df_dy,
                              Matrix<Real>& dg_dx, Matrix<Real>& dg_dy) const;
};
```

## Solution Access

```cpp
DAESolution& sol = result.solution;

// Time points
Real t = sol.getTValue(step);

// Differential variables
Real x_i = sol.getXValue(step, i);
Vector<Real> x = sol.getXValuesAtEnd();

// Algebraic variables
Real y_j = sol.getYValue(step, j);
Vector<Real> y = sol.getYValuesAtEnd();

// Interpolation
LinearInterpRealFunc xInterp = sol.getXAsLinInterp(0);
CubicSplineInterpRealFunc xSpline = sol.getXAsSplineInterp(0);
Real x_at_0_5 = xInterp(0.5);
```

## Example: Simple Linear DAE

```cpp
// System: dx/dt = -x + y, 0 = x + y - 1
// Analytical: x(t) = 0.5 + (x0 - 0.5)*exp(-2t), y = 1 - x

class SimpleLinearDAE : public IODESystemDAEWithJacobian {
public:
    int getDiffDim() const override { return 1; }
    int getAlgDim() const override { return 1; }

    void diffEqs(Real t, const Vector<Real>& x, const Vector<Real>& y,
                 Vector<Real>& dxdt) const override {
        dxdt[0] = -x[0] + y[0];
    }

    void algConstraints(Real t, const Vector<Real>& x, const Vector<Real>& y,
                        Vector<Real>& g) const override {
        g[0] = x[0] + y[0] - 1.0;
    }

    void jacobian_fx(Real t, const Vector<Real>& x, const Vector<Real>& y,
                     Matrix<Real>& df_dx) const override {
        df_dx(0, 0) = -1.0;
    }

    void jacobian_fy(Real t, const Vector<Real>& x, const Vector<Real>& y,
                     Matrix<Real>& df_dy) const override {
        df_dy(0, 0) = 1.0;
    }

    void jacobian_gx(Real t, const Vector<Real>& x, const Vector<Real>& y,
                     Matrix<Real>& dg_dx) const override {
        dg_dx(0, 0) = 1.0;
    }

    void jacobian_gy(Real t, const Vector<Real>& x, const Vector<Real>& y,
                     Matrix<Real>& dg_dy) const override {
        dg_dy(0, 0) = 1.0;  // Nonsingular → Index-1
    }
};

// Usage
SimpleLinearDAE system;
Vector<Real> x0(1), y0(1);
x0[0] = 0.8;
y0[0] = 0.2;  // Consistent: x + y = 1

DAESolverConfig config;
config.step_size = 0.01;

auto result = SolveDAEBDF2(system, 0.0, x0, y0, 1.0, config);

Real x_final = result.solution.getXValuesAtEnd()[0];
Real x_analytical = 0.5 + 0.3 * std::exp(-2.0);
// x_final ≈ 0.5406 (close to analytical solution)
```

## Mathematical Background

### Newton Iteration for Implicit DAE Step

At each time step, the implicit solver solves the coupled nonlinear system:

```
x_{n+1} - x_n - h*f(t_{n+1}, x_{n+1}, y_{n+1}) = 0
g(t_{n+1}, x_{n+1}, y_{n+1}) = 0
```

Using Newton iteration with the augmented Jacobian:

```
[I - h*∂f/∂x,  -h*∂f/∂y] [Δx]   [-F]
[  ∂g/∂x,        ∂g/∂y ] [Δy] = [-G]
```

Where:
- F = x_new - x_old - h*f (differential residual)
- G = g (algebraic residual)

### BDF2 Formula

For improved accuracy, BDF2 uses:

```
x_{n+1} = (4/3)*x_n - (1/3)*x_{n-1} + (2/3)*h*f(t_{n+1}, x_{n+1}, y_{n+1})
```

This provides second-order accuracy while maintaining A-stability.

## Index-1 Requirement

For the solver to work correctly, the algebraic Jacobian ∂g/∂y must be nonsingular. This ensures:

1. Unique solution for y given x at each time
2. Well-defined Newton iteration
3. Proper constraint satisfaction

If ∂g/∂y is singular, the system is higher-index and requires index reduction techniques (not currently supported).

## See Also

- [ODE Solvers](ODE_SOLVERS.md) - For ordinary differential equations
- [Linear Algebra](../base/LINEAR_ALGEBRA.md) - Matrix operations used internally
- [Numerical Analysis](../core/NUMERICAL_ANALYSIS.md) - Newton iteration and convergence

## References

1. Ascher, U.M., Petzold, L.R. (1998). *Computer Methods for Ordinary Differential Equations and Differential-Algebraic Equations*. SIAM.
2. Hairer, E., Wanner, G. (1996). *Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems*. Springer.
3. Brenan, K.E., Campbell, S.L., Petzold, L.R. (1996). *Numerical Solution of Initial-Value Problems in Differential-Algebraic Equations*. SIAM.
