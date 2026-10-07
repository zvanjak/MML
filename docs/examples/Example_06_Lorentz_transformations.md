# Example 06: Lorentz Transformations - Special Relativity

**Einstein's mind-bending physics, computed!** 🚀⏰

<p align="center">
  <img src="../images/examples/lorentz_twin_paradox.png" alt="Twin Paradox Worldlines" width="600">
</p>

## 📖 Overview

This self-contained example demonstrates **Special Relativity** while using MML for vectors, numerical differentiation, metric tensors, integration, interpolation, serialization, and visualization:

- **Time dilation** - Moving clocks tick slower
- **Length contraction** - Moving objects are shorter  
- **Twin Paradox** - The traveling twin ages less!
- **Spacetime diagrams** - Visualize worldlines
- **Accelerated journey** - Model acceleration, coasting, turnaround, and return
- **Proper-time verification** - Integrate the Minkowski metric numerically

> *"The distinction between past, present, and future is only a stubbornly persistent illusion."*  
> — Albert Einstein

**Source:** `src/examples/06_Lorentz_transformations/`

## 🎯 Scenarios

### Scenario 1: Time Dilation ⏰

**Moving clocks run SLOW.** A spaceship traveling at various speeds:

```
Relationship: τ = T / γ  where γ = 1/√(1 - v²/c²)

  Speed (v/c)    γ factor    Earth: 1 year    Ship: τ (years)
  -----------------------------------------------------------
     0.10        1.0050         1.0000          0.9950
     0.50        1.1547         1.0000          0.8660
     0.80        1.6667         1.0000          0.6000
     0.90        2.2942         1.0000          0.4359
     0.95        3.2026         1.0000          0.3122
     0.99        7.0888         1.0000          0.1411
     0.999      22.3663         1.0000          0.0447

  → At 99% of c, only 0.14 years pass on the ship!
  → At 99.9% of c, only 0.04 years pass!
```

This is NOT an illusion — it's a REAL physical effect, confirmed experimentally with atomic clocks on airplanes!

### Scenario 2: Length Contraction 📏

**Moving objects are SHORTER.** A 100-meter spaceship at various speeds:

```
Relationship: L = L₀ / γ = L₀ × √(1 - v²/c²)

  Speed (v/c)    γ factor    Measured Length (m)
  ------------------------------------------------
     0.00        1.00          100.00
     0.10        1.01           99.50
     0.50        1.15           86.60
     0.80        1.67           60.00
     0.90        2.29           43.59
     0.95        3.20           31.22
     0.99        7.09           14.11

  → At 90% of c, the 100m ship appears only 43.6m long!
  → At 99% of c, it's just 14.1m!
```

### Scenario 3: Twin Paradox 👬

The classic thought experiment: One twin stays on Earth, the other travels to Alpha Centauri and back.

```
Setup:
  - Twin A stays on Earth
  - Twin B travels to Alpha Centauri (4 light-years) and back
  - Ship speed: 0.8c (80% of light speed)
  - Lorentz factor γ: 1.67

Outbound Trip (Earth → Alpha Centauri):
  Earth time:  5.00 years
  Ship time:   3.00 years

Return Trip (Alpha Centauri → Earth):
  Earth time:  5.00 years
  Ship time:   3.00 years

═══════════════════════════════════════════════════════════════
  TOTAL JOURNEY:
    Twin A (Earth):  aged 10.00 years
    Twin B (Ship):   aged  6.00 years
    Difference:      4.00 years younger!
═══════════════════════════════════════════════════════════════
```

**Why is this NOT a paradox?** The situation is NOT symmetric! Twin B accelerates and decelerates, breaking the equivalence of inertial frames.

### Scenario 4: Spacetime Worldlines

Visualize the twins' paths through spacetime (Minkowski diagram):

- **Twin A:** Vertical line (stays at x=0, time passes)
- **Twin B:** V-shaped path (travels out and back)
- **Aging curves:** Proper time accumulated by each twin

### Scenario 5: Realistic Twin Paradox

The final simulation sends Twin B to a star 10 light-years away and back in seven phases:

1. Accelerate at approximately 1g from rest to 0.8c
2. Coast at 0.8c
3. Decelerate to rest near Star B
4. Complete a semicircular orbit around Star B
5. Accelerate to 0.8c toward Earth
6. Coast home
7. Decelerate to rest at Earth

The analytical model predicts approximately 30.11 years for Twin A and 20.06 years for Twin B, a difference of 10.05 years. It generates a 3D spacetime trajectory and an aging comparison plot.

As an independent check, the example represents Twin B's path as a 4D worldline, differentiates it numerically, and integrates proper time phase by phase:

$$\tau = \int \sqrt{-g_{\mu\nu}\frac{dx^\mu}{dt}\frac{dx^\nu}{dt}}\,dt$$

The numerical path integral is then compared with the analytical proper time for every phase and for the complete journey.

## ⚛️ Physics

### The Lorentz Transformation

Between two inertial frames moving at relative velocity v:

$$t' = \gamma(t - vx/c^2)$$
$$x' = \gamma(x - vt)$$
$$y' = y$$
$$z' = z$$

Where the **Lorentz factor:**

$$\gamma = \frac{1}{\sqrt{1 - v^2/c^2}}$$

### Natural Units

The example uses **natural units** where c = 1:
- Time measured in years
- Distance measured in light-years
- Velocity as fraction of c

This simplifies the equations:

$$\gamma = \frac{1}{\sqrt{1 - v^2}}$$

### Invariant Interval

The spacetime interval is invariant under Lorentz transformations:

$$ds^2 = c^2 dt^2 - dx^2 - dy^2 - dz^2$$

- **Timelike** (ds² > 0): Causal connection possible
- **Spacelike** (ds² < 0): No causal connection
- **Lightlike** (ds² = 0): Light path

## 🔧 MML Features Used

| Feature | Usage |
|---------|-------|
| `Vector4Minkowski` | 4-vector in Minkowski spacetime |
| `MetricTensorMinkowski` | Evaluate the spacetime interval along the worldline |
| `Derivation::DeriveCurve` | Numerically obtain the worldline tangent |
| `IntegrateTrap` | Independently integrate elapsed proper time |
| `LinearInterpRealFunc` | Smooth worldline interpolation |
| `Visualizer` | 2D aging plots and 3D spacetime trajectories |

The elementary x-axis Lorentz boost is deliberately implemented locally in `main.cpp`. It demonstrates how an MML user can combine a small domain-specific formula with MML's general numerical machinery without making special relativity part of the core library.

## 📁 Key Files

```
src/examples/06_Lorentz_transformations/
└── main.cpp              # Five demonstrations and local x-axis boost
```

## 🏃 Running

```bash
# Build
cmake --build build --config Release --target Example06_LorentzTransform

# Run
./build/src/examples/Release/Example06_LorentzTransform
```

## 📊 Sample Output

```
================================================================
     MML LORENTZ TRANSFORMATIONS - SPECIAL RELATIVITY DEMO
================================================================

   "The distinction between past, present, and future is only
    a stubbornly persistent illusion." - Albert Einstein

   Using natural units where c = 1 (speed of light).
   Time in years, distance in light-years.
================================================================

======================================================================
  SCENARIO 3: The Twin Paradox
======================================================================

Using the example's local LorentzBoostXAxis implementation:
  Event in Earth frame: t=10 years, x=0
  Same event in ship frame: t'=16.67 years, x'=-13.33

✓ Twin Paradox demonstration complete!
```

## 🌟 Key Insights

### Experimental Confirmation

These effects are **REAL**, not just theoretical:

1. **Hafele-Keating experiment (1971):** Atomic clocks flown around the world showed time dilation matching relativity
2. **GPS satellites:** Must correct for relativistic time dilation to maintain accuracy
3. **Muon experiments:** High-speed muons from cosmic rays live longer (from Earth's perspective) due to time dilation

### Common Misconceptions

| Misconception | Reality |
|---------------|---------|
| "It's just an illusion" | No, it's a real physical effect |
| "Only affects light" | Affects ALL objects and clocks |
| "Symmetric paradox" | Acceleration breaks symmetry |
| "Violates causality" | Preserves causality through light cone structure |

## 🎓 Learning Points

1. **Relativity of simultaneity:** Events simultaneous in one frame are NOT simultaneous in another
2. **Invariant speed of light:** c is the same for all observers — leads to all other effects
3. **Proper time:** The time measured by a clock traveling with an object (shortest path through spacetime)
4. **Composable tools:** A localized physics formula can build on MML's general vectors, tensors, differentiation, and integration

## 🔗 Related Topics

- **Example 02:** Double Pendulum (Edward Lorenz worked on atmospheric physics before discovering chaos)
- **General Relativity:** Curved spacetime, gravity as geometry (beyond this example)

## 📚 References

- Einstein, A. (1905). "On the Electrodynamics of Moving Bodies"
- Taylor & Wheeler, *Spacetime Physics* - Excellent introduction
- Landau & Lifshitz, *Classical Theory of Fields* - Advanced treatment
- Misner, Thorne & Wheeler, *Gravitation* - The bible of general relativity
