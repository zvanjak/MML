# Example 02: Double Pendulum - Chaos Theory in Action

**Where deterministic physics meets the unpredictable!**

<p align="center">
  <img src="../images/examples/double_pendulum_butterfly.png" alt="Butterfly Effect" width="600">
</p>

## 📖 Overview

The double pendulum is a **canonical example of deterministic chaos**. The equations are completely known, yet long-term prediction is impossible because infinitely small differences in initial conditions lead to completely different outcomes.

> *"Does the flap of a butterfly's wings in Brazil set off a tornado in Texas?"*  
> — Edward Lorenz, 1972

**Source:** `src/examples/02_double_pendulum/`

## 🎯 Scenarios

### Scenario 1: Chaotic Motion

Watch a single double pendulum evolve from a dramatic initial position (both arms at 135°).

```
Initial Conditions:
  θ₁ = 135° (upper arm)
  θ₂ = 135° (lower arm)
  ω₁ = ω₂ = 0 (released from rest)

Energy Conservation Check:
  Initial energy: -9.633814 J
  Final energy:   -9.633814 J
  Energy drift:   0.000001%   ← RK4 is EXCELLENT!
```

**Outputs:**
- `double_pendulum_angles.txt` - θ₁(t) and θ₂(t)
- `double_pendulum_phase_space.txt` - Chaotic attractor visualization
- `double_pendulum_trajectory.txt` - Full trajectory for animation

### Scenario 2: Butterfly Effect 🦋

**THE ESSENCE OF CHAOS:** Two pendulums starting with just **0.001°** difference!

```
Pendulum A: θ₁ = 90.000°, θ₂ = 90°
Pendulum B: θ₁ = 90.001°, θ₂ = 90°
                    └── JUST 0.001° !!

Divergence Over Time:
  Time (s)    Angle Difference
  --------------------------------
      0.0           0.0010 deg
      2.0           0.0089 deg
      5.0           0.4721 deg
      8.0          47.2314 deg   ← 47,000× amplification!!!
```

**The difference grows EXPONENTIALLY!** This is the butterfly effect in action.

### Scenario 3: Energy Conservation Check

Verify that numerical integration preserves the conserved quantities:

- Total mechanical energy E = T + V should remain constant
- RK4 typically achieves <0.001% drift over 10 seconds
- This validates our numerical approach

## ⚛️ Physics

### The System

A double pendulum consists of:
- **Mass 1** (m₁) at the end of rod of length l₁
- **Mass 2** (m₂) at the end of rod of length l₂, attached to m₁

### Equations of Motion

From Lagrangian mechanics (L = T - V):

**Angular accelerations:**

$$\ddot{\theta}_1 = \frac{-g(2m_1+m_2)\sin\theta_1 - m_2 g\sin(\theta_1-2\theta_2) - 2\sin(\theta_1-\theta_2)m_2[\dot{\theta}_2^2 l_2 + \dot{\theta}_1^2 l_1\cos(\theta_1-\theta_2)]}{l_1[2m_1+m_2-m_2\cos(2\theta_1-2\theta_2)]}$$

$$\ddot{\theta}_2 = \frac{2\sin(\theta_1-\theta_2)[\dot{\theta}_1^2 l_1(m_1+m_2) + g(m_1+m_2)\cos\theta_1 + \dot{\theta}_2^2 l_2 m_2\cos(\theta_1-\theta_2)]}{l_2[2m_1+m_2-m_2\cos(2\theta_1-2\theta_2)]}$$

These are **nonlinear coupled ODEs** — no closed-form solution exists!

### Total Energy (Conserved)

$$E = \frac{1}{2}m_1 l_1^2 \dot{\theta}_1^2 + \frac{1}{2}m_2[l_1^2\dot{\theta}_1^2 + l_2^2\dot{\theta}_2^2 + 2l_1 l_2\dot{\theta}_1\dot{\theta}_2\cos(\theta_1-\theta_2)]$$
$$- (m_1+m_2)gl_1\cos\theta_1 - m_2 g l_2\cos\theta_2$$

## 🔧 MML Features Used

| Feature | Usage |
|---------|-------|
| `IODESystem` | Define double pendulum dynamics |
| `ODESystemFixedStepSolver` | Fixed-step integration |
| `StepCalculators::RK4_Basic` | 4th-order Runge-Kutta |
| `PolynomInterpRealFunc` | Smooth interpolation for plotting |
| `SplineInterpParametricCurve` | Phase space trajectory |
| `Visualizer` | Multi-function plotting |

## 📁 Key Files

```
src/examples/02_double_pendulum/
├── main.cpp              # Demo scenarios
├── DoublePendulum.h      # Self-contained physics
└── CMakeLists.txt
```

## 🏃 Running

```bash
# Build
cmake --build build --target Example02_DoublePendulum

# Run
./build/src/examples/Release/Example02_DoublePendulum
```

## 📊 Sample Output

```
======================================================================
         MML DOUBLE PENDULUM - CHAOS DEMONSTRATION
======================================================================

   The double pendulum is a classic example of deterministic chaos.
   The motion is completely determined by the equations, yet
   impossible to predict long-term due to extreme sensitivity
   to initial conditions.

   "Does the flap of a butterfly's wings in Brazil set off
    a tornado in Texas?" - Edward Lorenz, 1972

======================================================================
  SCENARIO 1: Chaotic Double Pendulum Motion
======================================================================

Physical Parameters:
  Mass 1: 1 kg, Length 1: 1 m
  Mass 2: 1 kg, Length 2: 1 m

Initial Conditions:
  Theta 1: 135 degrees
  Theta 2: 135 degrees
  (Both arms pointing up-left - dramatic start!)

Solving for 10 seconds... Done! (2001 time steps)

Energy Conservation Check:
  Initial energy: -9.633814 J
  Final energy:   -9.633814 J
  Energy drift:   0.000001%

✓ Single trajectory complete!
```

## 🌀 What is Chaos?

**Chaos ≠ Randomness!**

| Property | Random | Chaotic |
|----------|--------|---------|
| Deterministic | No | **Yes** |
| Predictable (short-term) | No | **Yes** |
| Predictable (long-term) | No | No |
| Sensitive to initial conditions | N/A | **EXTREMELY** |

**Chaos is deterministic unpredictability.** The equations completely determine the future, but we can never know initial conditions with infinite precision.

## 🎓 Learning Points

1. **Chaos from simplicity:** Just 2 coupled oscillators → infinite complexity
2. **Exponential divergence:** Errors grow as e^(λt) where λ is the Lyapunov exponent
3. **Phase space:** The "attractor" shows structure in apparent randomness
4. **Energy conservation:** Validates numerical accuracy despite chaotic behavior
5. **MML workflow:** Same ODE framework handles simple and chaotic systems

## 🔗 Related Examples

- **Example 00:** N-Body Gravity (another chaotic system!)
- **Example 04:** Lorentz Transformations (Edward Lorenz discovered chaos!)

## 📚 References

- Strogatz, *Nonlinear Dynamics and Chaos* - THE textbook on chaos
- Lorenz (1963), "Deterministic Nonperiodic Flow" - Original chaos paper
- Wikipedia: "Double Pendulum" - Good mathematical treatment
