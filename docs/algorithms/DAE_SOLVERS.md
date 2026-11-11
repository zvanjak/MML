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
- **Index-1**: The algebraic Jacobian ∂g/∂y is nonsingular
- **Index-2+**: Higher-order DAEs requiring more derivatives to reduce

The MML solvers handle **index-1 systems** where ∂g/∂y is nonsingular.

The ordinary Cartesian pendulum position constraint `x*x + y*y - L*L = 0`
does not depend on the algebraic tension, so its algebraic Jacobian is zero.
That higher-index formulation requires index reduction; it is not an index-1
example merely because it implements the interface.

`IODESystemDAEWithJacobian` and `IODESystemDAEParametrized` share a virtual
`IODESystemDAE` base. A model can compose analytic Jacobians, parameters and
`IODESystemDAEWithEvents` without ambiguous ordinary-DAE conversions.

### Function-Pointer Adapter

`DAESystemWithJacobian` implements `IODESystemDAEWithJacobian` directly and can
be passed to both ordinary and Jacobian-capable DAE interface consumers,
including the implicit solvers below. Its differential/algebraic callbacks
reuse an internal `DAESystem`; there is one unambiguous interface chain.

**Compatibility note:** it is no longer a concrete `DAESystem` subclass.
Use `IODESystemDAE&` for ordinary DAE consumers instead of a concrete
`DAESystem&` upcast. Constructor arguments and evaluation methods are unchanged.

Configured construction requires both equation callbacks and all four Jacobian
callbacks; missing callbacks throw `ArgumentError`. Default construction remains
available for an empty adapter, but attempting any equation/Jacobian evaluation
throws `NotImplementedError` rather than silently leaving output buffers unchanged.

## Solver Methods

| Solver | Order | Stability | Step control | Typical use |
|---|---:|---|---|---|
| `SolveDAEBackwardEuler` | 1 | L-stable | Fixed | Robust baseline and very stiff systems |
| `SolveDAEBackwardEulerAdaptive` | 1 | L-stable | Adaptive step doubling | Variable dynamics and tolerance control |
| `SolveDAEBDF2` | 2 | A-stable | Fixed | General moderately stiff index-1 DAEs |
| `SolveDAEBDF4` | 4 | A(α)-stable | Fixed | Higher accuracy when the solution is smooth |
| `SolveDAELinearlyImplicit` | 2 | A-stable, not L-stable | Fixed | Linearly-implicit trapezoid (frozen Jacobian, no Newton); formerly misnamed `SolveDAERODAS` |
| `SolveDAERODAS4` | 4 | L-stable | Adaptive embedded error control | Strong stiff damping with a genuine six-stage Rosenbrock 4(3) pair |
| `SolveDAERadauIIA` | 5 | L-stable | Fixed | High-accuracy stiff integration |
| `SolveDAEBackwardEulerWithEvents` | 1 | L-stable | Fixed with event restarts | Hybrid systems and state discontinuities |

### Backward Euler (1st Order)

```cpp
#include <mml/algorithms/DAESolvers.h>
#include <mml/interfaces/IODESystemDAE.h>

// x' = -y, 0 = y - rate*x; g_y = 1, so this is index-1.
class ParametrizedDecay : public IODESystemDAEWithJacobian, public IODESystemDAEParametrized {
    Real rate = 1;
public:
    int getDiffDim() const override { return 1; }
    int getAlgDim() const override { return 1; }
    int getNumParam() const override { return 1; }
    Real getParam(int i) const override {
        if (i != 0) throw IndexError("Invalid parameter index");
        return rate;
    }
    void setParam(int i, Real value) override {
        if (i != 0) throw IndexError("Invalid parameter index");
        rate = value;
    }
    void diffEqs(Real, const Vector<Real>&, const Vector<Real>& y, Vector<Real>& f) const override { f[0] = -y[0]; }
    void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override { g[0] = y[0] - rate * x[0]; }
    void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 0; }
    void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = -1; }
    void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = -rate; }
    void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 1; }
};

ParametrizedDecay system;
system.setParams(Vector<Real>{2});
Vector<Real> x0{1}, y0{2}; // Consistent with y = rate*x.
DAESolverConfig config;
config.step_size = Real(1) / 8;
const auto result = SolveDAEBackwardEuler(system, 0, x0, y0, 1, config);
if (result.status != AlgorithmStatus::Success) throw NumericalMethodError(result.error_message);
```

The reduced exact solution is `x(t)=exp(-2*t)`, `y(t)=2*x(t)`. For eight
Backward Euler steps the independent discrete reference is `(1+2/8)^(-8)`.
`MML_DocsApp --interface-contracts` runs this example and verifies that recurrence,
the shared abstract base and `g_y=1`; see
[`docs_demo_interfaces.cpp`](../../src/docs_demos/core/docs_demo_interfaces.cpp).
The finite, nonnegative decay rate is an example input, not a new restriction on
all parametrized interfaces.

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

The shared implicit Newton solver and `ComputeConsistentICDetailed` check the
residual of the state after the final permitted correction. Their iteration
budgets limit corrections, not the final convergence check; diagnostic
residuals describe the returned state even when the budget is exhausted.
`iterations` counts attempted corrections. A zero budget in these direct
helpers checks the initial state without modifying it. The integration drivers
still require a positive `max_newton_iter` setting. Non-finite final states or
residuals remain explicit failures rather than being accepted at the budget
boundary.

### BDF4, linearly-implicit trapezoid, and Radau IIA

```cpp
auto bdf4 = SolveDAEBDF4(system, t0, x0, y0, tEnd, config);
auto linImpl = SolveDAELinearlyImplicit(system, t0, x0, y0, tEnd, config);
auto radau = SolveDAERadauIIA(system, t0, x0, y0, tEnd, config);
```

BDF4 builds its first three history intervals with the existing fifth-order,
L-stable Radau IIA solver, then uses BDF4 on the uniform grid. A materially
clipped final interval also uses Radau IIA, not a first-order Backward Euler
restart. This preserves fourth-order global accuracy on smooth index-1
problems, with Newton tolerances tight enough to resolve the discretization
error. A fixed number of BDF2 substeps starting with Backward Euler would
instead leave second-order starting error, regardless of how small the main
step becomes.

Each Radau history/endpoint interval uses a coupled three-stage system and is
more expensive than an ordinary BDF4 step. Its Newton/Jacobian/LU/solve work is
included in the BDF4 diagnostics, while `max_steps` and the saved trajectory
still count accepted main intervals, not stages. Failed candidates leave
history and `final_constraint_norm` at the last accepted state and increment
the rejection counters. Invalid used controls are rejected before equation
evaluation; adaptive-only controls do not affect this fixed-step method.
If the Radau starter cannot complete its requested interval, BDF4 reports
`Stalled`/`StepSizeUnderflow` rather than storing an unchanged candidate or
returning a success-shaped partial trajectory.

The linearly-implicit trapezoid uses a frozen augmented
Jacobian and linear predictor/corrector solves, so its Newton iteration count is
zero. The three-stage Radau IIA method solves a coupled stage system and is the
most expensive method per step, but provides fifth-order accuracy and strong
damping for stiff systems.

The linearly-implicit trapezoid is **second order** for smooth index-1 problems.
Its predictor is followed by one correction of the trapezoidal residual, using
the same frozen matrix. The correction subtracts the predictor increment and is
applied to the predicted state, not the old state. One Jacobian evaluation, one
LU factorization and two linear solves are used per successful step; there is
no Newton iteration loop.

On the linear test equation its amplification is
$R(z)=(1+z/2)/(1-z/2)$. This is **A-stable, not L-stable**:
$R(z)$ tends to $-1$, not zero, as $z$ tends to negative infinity. Unresolved
stiff modes can therefore alternate sign without useful damping. For strong
stiff damping use Backward Euler, Radau IIA or the genuine
[adaptive RODAS4 method](Differential_equations_solvers.md#genuine-rodas4-for-index-1-daes).
The deprecated `SolveDAERODAS` alias still selects the repaired trapezoid; it is
not redirected to RODAS4.

Nonlinear constraints are checked against `constraint_tol` after the single
correction, not projected by further iterations. A sufficiently small fixed
step may be necessary. Failed candidates do not enter the trajectory;
rejection counters are synchronized and `final_constraint_norm` retains the
last accepted state's norm. Only the fixed-step controls are required:
positive finite `step_size` and `constraint_tol`, a positive `max_steps` budget,
distinct finite endpoints and matching dimensions. Invalid controls throw
`ArgumentError` before equation evaluation. Newton/adaptive-only controls do
not change this method.

### Adaptive Backward Euler

```cpp
DAESolverConfig config;
config.step_size = 0.1;       // Initial step
config.abs_tolerance = 1e-8;
config.rel_tolerance = 1e-6;
config.min_step_size = 1e-8;
config.max_step_size = 0.5;

DAESolverResult result =
    SolveDAEBackwardEulerAdaptive(system, 0.0, x0, y0, 10.0, config);
```

The adaptive method compares one full Backward Euler step with two half steps,
rejects steps whose normalized error exceeds one, and adjusts the next step with
a bounded safety-factor controller. The accepted solution uses the two-half-step
result.

Excessive error is **never** accepted just because the minimum step has been
reached. An unattainable tolerance returns `AlgorithmStatus::NumericalInstability`
with `DAEFailureReason::StepSizeUnderflow` and a diagnostic. This also applies to
a final endpoint-clipped step shorter than `min_step_size`. Rejected candidates
do not enter the saved trajectory; accepted/rejected counters remain consistent
with its statistics, and `final_constraint_norm` describes the last accepted
state. `last_error_norm` and `max_error_norm` report scaled step-doubling errors,
including rejected estimates. Non-finite error estimates or overflowing scales
produce an explicit numerical failure rather than zero error.

Before solving, adaptive Backward Euler requires distinct finite endpoints,
matching state dimensions, positive finite step sizes and ordered bounds,
positive finite absolute/Newton/constraint tolerances, a nonnegative finite
relative tolerance, positive iteration/attempt budgets, and
`0 < step_safety < 1`. Invalid configuration throws `ArgumentError` before
step-size clamping or equation evaluation. Both accepted and rejected attempts
consume `max_steps`; budget exhaustion preserves the accepted partial solution.

## Configuration

```cpp
struct DAESolverConfig {
    Real step_size = 0.01;          // Fixed step or adaptive initial step
    int max_newton_iter = 20;       // Newton iterations per step
    Real newton_tol = 1e-8;         // Newton convergence tolerance
    Real constraint_tol = 1e-10;    // Constraint satisfaction tolerance
    int max_steps = 100000;         // Maximum integration steps
    Real abs_tolerance = 1e-10;     // Adaptive absolute tolerance
    Real rel_tolerance = 1e-8;      // Adaptive relative tolerance
    Real min_step_size = 1e-8;      // Adaptive lower step bound
    Real max_step_size = 1.0;       // Adaptive upper step bound
    Real step_safety = 0.9;         // Adaptive controller safety factor
};

// Preset configurations
DAESolverConfig highPrec = DAESolverConfig::HighPrecision();
DAESolverConfig fast = DAESolverConfig::Fast();
```

## Result Structure

Trajectory storage is bounded independently of `max_steps`. The initial condition
uses one column, so 100,000 accepted intervals require 100,001 points.
There is no fixed 100,000-column cap: capacity is derived from Matrix's unchanged
100-million-element budget for each differential/algebraic state matrix, the
time-vector element budget and integer index representability.
All DAE drivers, including event integration, retain the stored partial solution
and return `AlgorithmStatus::ResourceLimitExceeded` together with
`DAEFailureReason::ResourceLimitExceeded` when the next point cannot be reserved.
Saved-point counts and accepted-step counters do not advance for that point.
Growth is transactional across time, differential and algebraic arrays.
Fixed BDF4, Radau IIA and linearly-implicit budgets count accepted intervals,
not the initial point.

See [trajectory storage](Differential_equations_solvers.md#trajectory-storage-and-resource-limits)
for allocation limits, direct-container exceptions and the runnable demo.

```cpp
struct DAESolverResult {
    DAESolution solution;           // Solution trajectory
    
    std::string algorithm_name;     // Solver identifier
    AlgorithmStatus status;         // Success, NumericalInstability, etc.
    std::string error_message;
    double elapsed_time_ms;
    
    int total_steps;
    int newton_iterations;
    int jacobian_evaluations;
    int accepted_steps;
    int rejected_steps;
    int last_step_newton_iterations;
    int max_step_newton_iterations;
    Real final_residual_norm;
    Real max_residual_norm;
    Real final_constraint_norm;
    Real max_constraint_violation;
    DAEFailureReason failure_reason;
};
```

`failure_reason` distinguishes inconsistent initial conditions, Newton failure,
singular augmented Jacobians, likely higher-index systems, maximum-step exits,
constraint violations, and non-finite states.

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

For diagnostics and damped Newton recovery, use the detailed API:

```cpp
ConsistentICResult ic = ComputeConsistentICDetailed(system, t0, x0, y0);
if (!ic.converged) {
    std::cerr << ic.message << " (residual " << ic.final_residual_norm << ")\n";
}
```

The detailed result reports iterations, line-search reductions, initial/final
residual norms, and a structured failure reason. A singular `dg/dy` is reported
as a likely higher-index or incorrectly formulated system.

## Numerical Jacobians

Analytical Jacobians are recommended when available. A system implementing only
`IODESystemDAE` can use all implicit solvers through the fourth-order finite-
difference adapter:

```cpp
MyDAEWithoutJacobians baseSystem;
DAESystemNumericalJacobian system(baseSystem);
auto result = SolveDAEBDF2(system, t0, x0, y0, tEnd, config);
```

The adapter computes `df/dx`, `df/dy`, `dg/dx`, and `dg/dy`, and supports the
combined `allJacobians` path. Automatic-differentiation and sparse Jacobian
policies are future scalability work; they are not required for current dense
index-1 solvers.

## Events And Reinitialization

Event-capable systems derive from `IODESystemDAEWithEvents`. Event functions can
depend on both differential and algebraic state and use the same
`EventDirection` and `EventAction` conventions as ODE events.

```cpp
class SwitchedDAE : public IODESystemDAEWithEvents {
    // Implement the DAE equations and four Jacobian blocks as usual.
public:
    int getNumEvents() const override { return 1; }

    Real eventFunction(int, Real t, const Vector<Real>& x,
                       const Vector<Real>& y) const override {
        return x[0] - 1.0;
    }

    EventDirection getEventDirection(int) const override {
        return EventDirection::Increasing;
    }

    EventAction getEventAction(int) const override {
        return EventAction::Restart;
    }

    void handleEvent(int, Real, Vector<Real>& x, Vector<Real>& y) const override {
        x[1] = -x[1];       // Apply the discrete state change.
        y[0] = 0.0;         // Initial guess; consistency is recomputed next.
    }
};

DAEEventConfig events;
events.event_tolerance = 1e-10;
auto result = SolveDAEBackwardEulerWithEvents(
    system, t0, x0, y0, tEnd, config, events);
```

Actions are:

- `Continue`: record the root and continue from the projected event state.
- `Stop`: record the root and terminate at the event time.
- `Restart`: invoke `handleEvent`, recompute consistent algebraic variables with
  damped Newton, validate constraints, and restart integration at the root.

`DAEEventResult` contains the integration diagnostics, chronological event
records with both `x` and `y`, terminal-event status, and final state. Event
states are projected onto the constraint manifold before use.

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

## Validation Guidance

For a new DAE model, validate more than the final differential state:

1. Verify or compute consistent initial conditions before integration.
2. Check `status`, `failure_reason`, and `error_message` before consuming output.
3. Monitor `max_constraint_violation` and `final_constraint_norm`.
4. Repeat with half the fixed step, or tighter adaptive tolerances, and compare.
5. For stiff models, compare at least two L-stable methods such as Backward Euler
    and Radau IIA.
6. Validate event direction, root time, post-event state, and constraints after
    every restart.

Core tests cover linear and nonlinear index-1 systems, stiff systems, electrical
circuits, constraint drift, all solver families, numerical Jacobians, adaptive
step rejection, singular/higher-index diagnostics, non-finite states, and DAE
event continuation/termination/reinitialization.

## See Also

- [ODE Solvers](ODE_SOLVERS.md) - For ordinary differential equations
- [Linear Algebra](../base/LINEAR_ALGEBRA.md) - Matrix operations used internally
- [Numerical Analysis](../core/NUMERICAL_ANALYSIS.md) - Newton iteration and convergence

## References

1. Ascher, U.M., Petzold, L.R. (1998). *Computer Methods for Ordinary Differential Equations and Differential-Algebraic Equations*. SIAM.
2. Hairer, E., Wanner, G. (1996). *Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems*. Springer.
3. Brenan, K.E., Campbell, S.L., Petzold, L.R. (1996). *Numerical Solution of Initial-Value Problems in Differential-Algebraic Equations*. SIAM.
