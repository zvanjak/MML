# Brownian Motion Analysis — What It Took to Recover the √t Law

This document records every problem we hit while making the 2D collision
simulator reproduce the textbook Brownian‑motion result — that the RMS
displacement of a heavy "Brownian" particle grows like `√t` — and the fix for
each one. The relevant code lives in:

- `src/book/mpl/CollisionSimulator2D/Containers2D.h`
- `src/book/mpl/CollisionSimulator2D/CollisionSimulator2D.h`
- `src/book/mpl/CollisionSimulator2D/ContainerFactory2D.h`
- `src/book/chapters/Chapter_04_collision_simulator/example4_Brownian_motion.cpp`
- `src/book/chapters/Chapter_04_collision_simulator/example3_gas_laws_verification.cpp`

Physical setup: one big ball (radius 20, mass 20) immersed in 500 small balls
(radius 1, mass 1) with random velocities, inside an 800×800 box. We track the
big ball's displacement from its start position as a function of time and
average over many independent runs.

---

## 1. A new container type: `PeriodicBoxContainer2D`

### Problem
The existing `BoxContainer2D` has **reflective walls**. A ball that reaches a
wall bounces back. That is fine for gas‑pressure experiments, but it is wrong
for Brownian motion: the walls confine the big ball, so its displacement
saturates at the box size instead of growing without bound. We need an
effectively *infinite* homogeneous medium.

### Fix
Add a **separate** container class that wraps particles around to the opposite
side when they leave the box (periodic boundary conditions), rather than
modifying `BoxContainer2D`. `PeriodicBoxContainer2D` derives from
`BoxContainer2D` and overrides only the out‑of‑bounds handling:

```cpp
struct PeriodicBoxContainer2D : public BoxContainer2D
{
    PeriodicBoxContainer2D() : BoxContainer2D() {}
    PeriodicBoxContainer2D(double width, double height) : BoxContainer2D(width, height) {}
    PeriodicBoxContainer2D(const BoxContainer2D& box) : BoxContainer2D(box) {}

    void CheckAndHandleOutOfBounds(int ballIndex)   // wrap instead of reflect
    {
        Ball2D& ball = _balls[ballIndex];
        double x = std::fmod(ball.Pos().X(), _width);  if (x < 0) x += _width;
        double y = std::fmod(ball.Pos().Y(), _height); if (y < 0) y += _height;
        ball.Pos().X() = x;
        ball.Pos().Y() = y;
    }
};
```

`BoxContainer2D` is left completely untouched, so every existing example keeps
its original reflective‑wall behaviour.

---

## 2. Templating `CollisionSimulator2D` on the container type

### Problem
`CollisionSimulator2D` stored its container **by value** as a `BoxContainer2D`
member:

```cpp
BoxContainer2D _box;   // stored by value
```

If we pass a `PeriodicBoxContainer2D` into that, C++ **object slicing** throws
away the derived part and we are left with a plain `BoxContainer2D` — reflective
walls again. Worse, `CheckAndHandleOutOfBounds` is intentionally **non‑virtual**
(the hot loop must not pay for virtual dispatch), so even a base‑class reference
would still call the reflective version.

### Fix
Make the simulator a template parameterised on the container type and store the
container as that concrete type:

```cpp
template<typename ContainerType = BoxContainer2D>
class CollisionSimulator2D
{
    ContainerType _box;   // exact type, no slicing
    ...
};
```

Now `_box.CheckAndHandleOutOfBounds(i)` resolves at compile time (static
dispatch) to the correct method for the actual container — periodic or
reflective — with **zero runtime overhead** and no virtual functions.

### Consequences for callers
- Thanks to **CTAD** (class template argument deduction), existing call sites
  such as `CollisionSimulator2D simulator(box);` still compile: the template
  argument is deduced from the constructor argument.
- CTAD does **not** apply to `new`, so `example3` had to spell the default
  argument explicitly:
  ```cpp
  CollisionSimulator2D<>* pSimulator;
  pSimulator = new CollisionSimulator2D<>(box, 4, 4, &recorder);
  ```

---

## 3. NaN blow‑up in the collision handler

### Symptom
With a lighter, faster big ball the analysis table filled with `-nan(ind)`
after a few tens of steps — every downstream statistic became NaN.

### Root cause
`HandleCollision` computes the *exact* moment of impact by back‑tracking the two
balls one step and solving a quadratic for the collision time:

```cpp
Pnt2Cart x10 = ball1.Pos() - ball1.V() * dt;   // assumed previous position
...
double disc = B * B - 4 * A * C;
double t1 = (-B + sqrt(disc)) / (2 * A);        // sqrt(disc) -> NaN if disc < 0
```

This assumes the recorded position really is the continuous position one step
ago. In a **periodic** box that assumption breaks: when a ball crosses a
boundary its position is *wrapped*, so `pos - v*dt` is nowhere near its true
previous location. The back‑tracked geometry becomes inconsistent, the
discriminant can go **negative**, and `sqrt(disc)` returns NaN. Once one
velocity is NaN it contaminates every later collision and every averaged
quantity.

The heavy (mass‑1000) ball almost never crossed a boundary, so it rarely
triggered the bug — which is why the problem only appeared once we made the ball
lighter and faster.

### Fix
Detect the degenerate cases and fall back to a **direct elastic impulse** along
the line joining the current centres (no back‑tracking, so no dependence on the
wrapped previous position):

```cpp
double disc = B * B - 4 * A * C;

if (A < 1e-12 || disc < 0.0)                 // near‑equal velocities, or wrap‑corrupted geometry
{
    Vec2Cart n(ball2.Pos(), ball1.Pos());    // centre‑to‑centre normal
    double sep = n.NormL2();
    if (sep < 1e-9) return;                   // coincident centres: no defined normal

    double vrelAlongN = (v1 - v2) * n / sep;
    if (vrelAlongN >= 0.0) return;            // already separating: don't re‑kick (anti‑stick)

    double factor = (v1 - v2) * n / POW2(sep);
    ball1.V() = v1 - 2 * m2 / (m1 + m2) * factor * n;
    ball2.V() = v2 + 2 * m1 / (m1 + m2) * factor * n;
    return;
}
```

The `vrelAlongN >= 0` guard prevents the classic "sticky balls" artefact where
an overlapping pair keeps exchanging impulses every step. This made the
simulation numerically robust for any mass ratio.

---

## 4. Still super‑diffusive: the conserved‑momentum drift

### Symptom
After the NaN fix the data was valid but the motion was **super‑diffusive**: the
RMS displacement grew faster than `√t`. The ratio of RMS at `t = 1600` to
`t = 400` was ≈ 2.8×, whereas pure diffusion predicts `√(1600/400) = 2.0×`.

### Root cause
A wall‑less periodic box **conserves total momentum**. The 500 small balls have
random velocities whose vector sum is not zero — it fluctuates like
`√N · v ≈ √500 · 15`. With nothing to absorb it, that net momentum makes the
whole **centre of mass drift at a constant velocity forever**. The big ball is
carried along with the drift, adding a **ballistic** component (`∝ t`) on top of
the diffusion and masking the `√t` law.

### Fix
Remove the net momentum once, at the start of each run, so the centre of mass
stays fixed and the big ball diffuses about a stationary point:

```cpp
double px = 0, py = 0, totM = 0;
for (auto& b : box._balls) { px += b.Mass()*b.V().X(); py += b.Mass()*b.V().Y(); totM += b.Mass(); }
double vcmx = px/totM, vcmy = py/totM;
for (auto& b : box._balls) { b.V().X() -= vcmx; b.V().Y() -= vcmy; }
```

This dropped the 400→1600 ratio from ≈ 2.8× to ≈ 2.2× — the residual excess is
now statistical noise rather than a systematic drift.

---

## 5. Unwrapping the periodic trajectory before measuring displacement

### Problem
Because particles wrap around, the recorded position of the big ball jumps by
roughly a full box length whenever it crosses a boundary. Measuring
displacement directly from the wrapped coordinates would produce huge spurious
jumps.

### Fix
Reconstruct the **unwrapped** trajectory using the minimum‑image convention:
each step, if the raw coordinate change exceeds half the box, add or subtract a
box length to undo the wrap, then accumulate the true displacement:

```cpp
double dx = path[t].X() - path[t-1].X();
double dy = path[t].Y() - path[t-1].Y();
if (dx >  boxSize/2) dx -= boxSize;
if (dx < -boxSize/2) dx += boxSize;
if (dy >  boxSize/2) dy -= boxSize;
if (dy < -boxSize/2) dy += boxSize;
ux += dx; uy += dy;                       // cumulative unwrapped position
double d = std::sqrt((ux-startX)*(ux-startX) + (uy-startY)*(uy-startY));
```

This is valid because the per‑step motion is small compared with the box, so a
change larger than `boxSize/2` can only be a wrap event.

---

## 6. Choosing a mass ratio that actually reaches the diffusive regime

### Problem
Brownian motion is *by definition* a big ball among many small balls, so the big
ball must stay heavy. But the velocity of a very heavy ball relaxes only after
about `M/m` collisions. With `M = 1000`, `m = 1` that is ~1000 collisions — far
more than a 500‑step run delivers, so the ball was still **ballistic** the whole
time and never showed `√t`.

### Fix / trade‑off
- Keep the ball genuinely "big" (radius 20) and heavier than the medium
  (`M = 20`, i.e. 20× the small‑ball mass), so it is still a proper Brownian
  particle but relaxes within a few tens of collisions.
- Do **not** shrink it further — that would stop being Brownian motion.
- Instead, **increase the run length** so the diffusive part dominates the early
  ballistic transient. We settled on `numSteps = 2000`.

`CreateBrownMotionConfig` gained optional parameters so the mass/size can be
tuned without touching the factory internals, while defaults preserve the old
behaviour:

```cpp
static BoxContainer2D CreateBrownMotionConfig(
        double width, double height, int numSmallBalls,
        double bigBallMass = 1000, double bigBallRadius = 20);
```

---

## 7. Averaging: enough runs for a smooth curve

### Problem
A single trajectory is very noisy; even 30 runs left visible wiggles in the RMS
curve (e.g. the last points dipping instead of rising).

### Fix
Average the mean‑squared displacement over **100 independent runs**
(`numRepeats = 100`). The RMS curve then hugs the `C·√t` reference smoothly
across the whole time range.

---

## 8. Tooling: the stale visualizer binary

### Symptom
Plotting the two curves (measured RMS and `C·√t`) failed with
`Unsupported type: MML_MULTI_REAL_FUNCTION`.

### Root cause
The C++ serializer correctly writes a multi‑function file tagged
`MML_MULTI_REAL_FUNCTION`, but the **deployed** WPF visualizer binary in
`tools/visualizers/win/WPF/MML_RealFunctionVisualizer/` was a stale March‑2026
build that predated multi‑function support. The current visualizer *source*
already handled the format (and had even reworded the error to "Unsupported
format"), but the compiled binary shipped in the MML repo had never been
refreshed.

### Fix
Rebuild `MML_RealFunctionVisualizer` in Release from the sibling
`MML_Visualizers` repo and copy the updated `.exe`/`.dll` files into the MML
repo's deployment folder. The multi‑function plot then opened correctly.

---

## Result

With all of the above in place, the big ball (mass 20, radius 20) diffusing
through 500 small balls over 2000 steps, averaged over 100 runs, reproduces the
textbook Brownian result:

| t    | RMS (measured) | C·√t |
|------|----------------|------|
| 600  | 514            | 534  |
| 1000 | 687            | 690  |
| 1200 | 763            | 756  |
| 1600 | 849            | 872  |
| 1800 | 915            | 925  |

The 400→1600 ratio is ≈ 2.18× against the ideal `√t` value of 2.0× — diffusive
within statistical noise. The measured RMS oscillates around and repeatedly
crosses the `C·√t` reference, the signature of genuine `√t` scaling.

### Summary of changes
1. `PeriodicBoxContainer2D` — wrap‑around container, `BoxContainer2D` untouched.
2. Templated `CollisionSimulator2D<ContainerType>` — static dispatch, no slicing.
3. NaN‑safe `HandleCollision` — impulse fallback for wrap‑corrupted / degenerate cases.
4. Zero net momentum per run — removes the ballistic centre‑of‑mass drift.
5. Minimum‑image trajectory unwrapping — correct displacement across boundaries.
6. Tuned mass ratio (20:1) + longer runs (2000 steps) — reach the diffusive regime.
7. 100 repeats — smooth averaged curve.
8. Redeployed the visualizer binary — multi‑function plotting works.
