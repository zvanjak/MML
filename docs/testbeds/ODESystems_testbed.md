# ODE Systems Test Bed Documentation

## Overview

The ODE Systems Test Beds provide comprehensive test infrastructure for validating MML's ordinary differential equation solvers. The framework includes two complementary test beds:

1. **ODESystemTestBed** - General ODE systems (linear, oscillatory, nonlinear, chaotic, physical)
2. **Stiff ODE Test Bed** - Classic stiff benchmark problems for implicit solvers

**Key Features:**
- Systems with **analytical solutions** for accuracy verification
- Systems with **end-point solutions** for integration testing  
- Systems **with Jacobians** for implicit solver testing
- Category-based organization and filtering
- Stiffness ratio metadata for solver selection

---

## Part 1: General ODE Systems Test Bed

### Core Files

| File | Purpose |
|------|---------|
| [diff_eq_systems_test_bed.h](../../test_data/diff_eq_systems_test_bed.h) | `ODESystemTestBed` class, interfaces, accessors |
| [diff_eq_systems_defs.h](../../test_data/diff_eq_systems_defs.h) | ODE system class definitions |

### Interfaces

```cpp
namespace MML::TestBeds {

/// Interface for ODE systems with known analytical solutions
/// Allows verification at ANY time point within the integration domain
class ITestODESystemWithSolution {
    virtual const IODESystem* getODESystem() const = 0;
    virtual Vector<Real> getInitCondLow() const = 0;
    virtual Vector<Real> getInitCondHigh() const = 0;
    virtual Vector<Real> getSolution(const Vector<Real>& initCond, Real t) const = 0;
};

/// Interface for ODE systems with known end-point solutions
/// Perfect for unit testing: integrate from t0 to tEnd and compare
class ITestODESystemWithEndSolution {
    virtual const IODESystem* getODESystem() const = 0;
    virtual Vector<Real> getInitialConditions() const = 0;
    virtual Real getStartTime() const = 0;
    virtual Real getEndTime() const = 0;
    virtual Vector<Real> getEndSolution() const = 0;
};

}
```

### System Categories

```cpp
enum class ODECategory {
    Linear,           // Linear systems (exact analytical solutions)
    Oscillatory,      // Harmonic oscillators (simple, damped)
    Nonlinear,        // Nonlinear oscillators, pendulum
    Population,       // Population dynamics (Lotka-Volterra, logistic)
    Chaotic,          // Lorenz, Rossler
    Physical,         // Projectile, Kepler
    SpecialFunction,  // Legendre, Bessel, Hermite, Laguerre
    Stiff             // Stiff systems (single example)
};
```

### Test Systems Catalog

#### 1. Linear Systems (with Analytical Solutions)

| # | Name | Equations | Dim | Solution |
|---|------|-----------|-----|----------|
| 1 | ExponentialDecayODE | y' = -ky | 1 | y(t) = y₀·e^(-kt) |
| 2 | ExponentialGrowthODE | y' = ky | 1 | y(t) = y₀·e^(kt) |
| 3 | Linear2DDistinctEigenODE | y₁' = -y₁+y₂, y₂' = -4y₂ | 2 | λ₁=-1, λ₂=-4 |
| 4 | Linear3DSystemODE | 3×3 linear system | 3 | λ = {1, 2, 2} |

#### 2. Oscillatory Systems (with Analytical Solutions)

| # | Name | Equations | Parameters | Notes |
|---|------|-----------|------------|-------|
| 5 | SimpleHarmonicOscillatorODE | x'' + ω²x = 0 | ω=1.0, 2.0 | Energy conserving |
| 6 | DampedHarmonicOscillatorODE | x'' + 2ζωx' + ω²x = 0 | ω=1.0, ζ=0.1 | Underdamped case |

#### 3. Nonlinear Oscillators (No General Analytical Solution)

| # | Name | Equations | Parameters | Behavior |
|---|------|-----------|------------|----------|
| 7 | VanDerPolODE | x'' - μ(1-x²)x' + x = 0 | μ=0.3 | Limit cycle |
| 8 | DuffingODE | x'' + δx' + αx + βx³ = γcos(ωt) | δ=0.3, α=1, β=0.2 | Nonlinear spring |
| 9 | SimplePendulumODE | θ'' + (g/L)sin(θ) = 0 | g/L=9.81 | Large angle |

#### 4. Population Dynamics

| # | Name | Equations | Parameters | Conserved Quantity |
|---|------|-----------|------------|-------------------|
| 10 | LotkaVolterraODE | y₁' = αy₁ - βy₁y₂, y₂' = δy₁y₂ - γy₂ | α=1.1, β=0.4, δ=0.1, γ=0.4 | First integral |
| 11 | LogisticGrowthODE | y' = ry(1 - y/K) | r=1.0, K=10.0 | Analytical solution |

#### 5. Chaotic Systems

| # | Name | Equations | Parameters | Notes |
|---|------|-----------|------------|-------|
| 12 | LorenzSystemODE | x'=σ(y-x), y'=x(ρ-z)-y, z'=xy-βz | σ=10, ρ=28, β=8/3 | Classic strange attractor |
| 13 | RosslerSystemODE | x'=-y-z, y'=x+ay, z'=b+z(x-c) | a=0.2, b=0.2, c=5.7 | Simpler chaotic system |

#### 6. Physical Systems

| # | Name | Equations | Dim | Conservation Laws |
|---|------|-----------|-----|-------------------|
| 14 | ProjectileWithDragODE | x'' = -(c/m)vₓ|v|, y'' = -g-(c/m)vᵧ|v| | 4 | None (dissipative) |
| 15 | KeplerProblemODE | r'' - rθ'² = -μ/r² | 4 | Energy, angular momentum |

#### 7. Special Function ODEs

| # | Name | Differential Equation | Order |
|---|------|----------------------|-------|
| 16 | LegendreODE | (1-x²)y'' - 2xy' + n(n+1)y = 0 | n=2 |
| 17 | LaguerreODE | xy'' + (1-x)y' + ny = 0 | n=3 |
| 18 | HermiteODE | y'' - 2xy' + 2ny = 0 | n=2 |
| 19 | BesselODE | x²y'' + xy' + (x² - n²)y = 0 | n=0 |

### Systems with End-Point Solutions

Pre-computed solutions for quick integration tests:

| Name | Category | IC | t₀ | tₑₙd | Description |
|------|----------|-----|-----|------|-------------|
| ExpDecay_t5 | Linear | [1.0] | 0 | 5 | Exponential decay |
| ExpGrowth_t2 | Linear | [1.0] | 0 | 2 | Exponential growth |
| Linear2D_t3 | Linear | [1.0, 3.0] | 0 | 3 | 2D linear system |
| Linear3D_t1 | Linear | [1,1,1] | 0 | 1 | 3D linear system |
| SHO_1period | Oscillatory | [1.0, 0.0] | 0 | 2π | One period |
| SHO_5periods | Oscillatory | [1.0, 0.0] | 0 | 10π | Energy test |
| DampedSHO_t10 | Oscillatory | [1.0, 0.0] | 0 | 10 | Damped oscillation |
| DampedSHO_t20 | Oscillatory | [1.0, 0.0] | 0 | 20 | Long-time decay |
| Logistic_t5 | Population | [1.0] | 0 | 5 | Logistic growth |
| Logistic_t10 | Population | [0.5] | 0 | 10 | Approach carrying capacity |

---

## Part 2: Stiff ODE Systems Test Bed

### Core Files

| File | Purpose |
|------|---------|
| [stiff_ode_test_bed.h](../../test_data/stiff_ode_test_bed.h) | `TestStiffODE` struct, retrieval functions |
| [stiff_ode_defs.h](../../test_data/stiff_ode_defs.h) | Stiff ODE system class definitions with Jacobians |

### Data Structure

```cpp
namespace MML::TestBeds {

/// Test case for stiff ODE solvers
struct TestStiffODE {
    std::string name;                          ///< Descriptive name
    int dimension;                             ///< System dimension
    
    // Factory functions
    std::function<std::unique_ptr<IODESystem>()> createSystem;
    std::function<std::unique_ptr<IODESystemWithJacobian>()> createStiffSystem;
    
    // Problem setup
    Vector<Real> initialCondition;             ///< y(0)
    Real tStart = 0.0;                         ///< Integration start
    Real tEnd = 1.0;                           ///< Integration end
    
    // Reference solution (if known)
    bool hasExactSolution = false;
    std::function<Vector<Real>(Real)> exactSolution;
    
    // Stiffness characteristics
    Real stiffnessRatio = 1.0;                 ///< Fastest/slowest eigenvalue
    std::string stiffnessCategory;             ///< mild, moderate, severe, extreme
    
    // Metadata
    std::string description;
    std::string reference;                     ///< Literature reference
    
    // Testing hints
    Real suggestedStepSize = 0.01;
    Real expectedAccuracy = 1e-6;
    int difficulty = 2;                        ///< 1=easy, 4=extreme
};

}
```

### Stiff Systems Catalog (9 total)

#### Chemical Kinetics

| # | Name | Dim | Stiffness | Ratio | Description |
|---|------|-----|-----------|-------|-------------|
| 1 | **Robertson** | 3 | Extreme | 10⁸ | Classic chemical kinetics A→B→C |
| 2 | **RobertsonScaled** | 2 | Extreme | 10⁸ | Conserved form (2D) |
| 3 | **Brusselator** | 2 | Moderate | ~5 | Autocatalytic oscillator |
| 4 | **Oregonator** | 3 | Severe | 10⁷ | Belousov-Zhabotinsky reaction |
| 5 | **E5_Decay** | 4 | Extreme | 6×10⁹ | Linear decay chain |

#### Physics/Biology

| # | Name | Dim | Stiffness | Description |
|---|------|-----|-----------|-------------|
| 6 | **VanDerPolStiff** | 2 | Severe | μ=1000, relaxation oscillator |
| 7 | **HIRES** | 8 | Moderate | Plant photomorphogenesis |
| 8 | **Pollution** | 10 | Moderate | Atmospheric chemistry |

#### With Exact Solution

| # | Name | Dim | Stiffness | Solution |
|---|------|-----|-----------|----------|
| 9 | **LinearStiff** | 3 | Severe (10⁴) | y(t) = [e⁻ᵗ, e⁻¹⁰⁰⁰ᵗ, e⁻¹⁰⁰⁰⁰ᵗ] |

### Stiffness Categories

| Category | Stiffness Ratio | Typical Step Restriction | Example |
|----------|----------------|-------------------------|---------|
| **Mild** | < 10³ | h < 0.1 | - |
| **Moderate** | 10³ - 10⁵ | h < 10⁻³ | Brusselator, HIRES |
| **Severe** | 10⁵ - 10⁸ | h < 10⁻⁵ | Oregonator, VanDerPol(μ=1000) |
| **Extreme** | > 10⁸ | h < 10⁻⁸ | Robertson, E5 |

---

## API Reference

### ODESystemTestBed Accessors

```cpp
class ODESystemTestBed {
public:
    // Counts
    static int numSystemsWithSolution();        // Systems with analytical solutions
    static int numSystemsWithEndSolution();     // Systems with end-point solutions
    static int numGeneralSystems();             // Systems without solutions (qualitative)
    static int numSystemsWithJacobian();        // Systems for implicit methods
    
    // Access by index
    static const IODESystem* getSystemWithSolution(int index);
    static Vector<Real> getSolution(int index, const Vector<Real>& ic, Real t);
    static const ODESystemWithEndSolutionEntry& getSystemWithEndSolution(int index);
    static const ODESystemEntry& getGeneralSystem(int index);
    static const IODESystemWithJacobian* getSystemWithJacobian(int index);
    
    // Access by name
    static const IODESystem* getSystemWithSolutionByName(const std::string& name);
    static const ODESystemWithEndSolutionEntry* getSystemWithEndSolutionByName(const std::string& name);
    static const IODESystem* getGeneralSystemByName(const std::string& name);
    
    // Access by category
    static std::vector<const ODESystemWithSolutionEntry*> getSystemsWithSolutionByCategory(ODECategory cat);
    static std::vector<const ODESystemWithEndSolutionEntry*> getSystemsWithEndSolutionByCategory(ODECategory cat);
    static std::vector<const ODESystemEntry*> getGeneralSystemsByCategory(ODECategory cat);
    
    // Filter by dimension
    static std::vector<const ODESystemWithEndSolutionEntry*> getSystemsWithEndSolutionByDimension(int dim);
    static std::vector<const ODESystemEntry*> getGeneralSystemsByDimension(int dim);
    
    // Convenience methods
    static std::vector<const ODESystemWithEndSolutionEntry*> getLinearSystems();
    static std::vector<const ODESystemWithEndSolutionEntry*> getOscillatorySystems();
    static std::vector<const ODESystemEntry*> getChaoticSystems();
    static std::vector<const ODESystemWithEndSolutionEntry*> get1DSystems();
    static std::vector<const ODESystemWithEndSolutionEntry*> get2DSystems();
};
```

### Stiff ODE Test Bed Functions

```cpp
namespace MML::TestBeds {

// Individual test case generators
TestStiffODE getRobertsonTest();
TestStiffODE getVanDerPolStiffTest(Real mu = 1000.0);
TestStiffODE getBrusselatorTest(Real A = 1.0, Real B = 3.0);
TestStiffODE getOregonatorTest();
TestStiffODE getHIRESTest();
TestStiffODE getRobertsonScaledTest();
TestStiffODE getE5Test();
TestStiffODE getPollutionTest();
TestStiffODE getLinearStiffTest();

// Collection getters
std::vector<TestStiffODE> getAllStiffODETests();           // All 9 tests
std::vector<TestStiffODE> getModerateStiffTests();         // Brusselator, HIRES, Linear
std::vector<TestStiffODE> getSevereStiffTests();           // Robertson, VanDerPol, Oregonator, E5
std::vector<TestStiffODE> getExactSolutionStiffTests();    // LinearStiff only
std::vector<TestStiffODE> getChemicalKineticsTests();      // Robertson, Brusselator, Oregonator, E5
std::vector<TestStiffODE> getSmallDimensionStiffTests();   // 2-3D systems
std::vector<TestStiffODE> getLargeDimensionStiffTests();   // HIRES (8D), Pollution (10D)

// Verification utilities
Real computeStiffODEError(const Vector<Real>& computed, const TestStiffODE& test, Real t);
Real checkConservation(const Vector<Real>& y, const TestStiffODE& test);
bool verifyNonNegative(const Vector<Real>& y);

}
```

---

## Usage Examples

### Example 1: Test Adaptive Solver with Analytical Solution

```cpp
#include "test_data/diff_eq_systems_test_bed.h"
#include "algorithms/ODEAdaptiveIntegrator.h"

using namespace MML;
using namespace MML::TestBeds;

void testAdaptiveSolver() {
    // Get Simple Harmonic Oscillator system
    auto& entry = *ODESystemTestBed::getSystemWithEndSolutionByName("SHO_1period");
    const IODESystem* sys = entry.system;
    
    // Integrate
    CashKarpIntegrator integrator(*sys);
    auto solution = integrator.integrate(
        entry.ic,           // Initial conditions [1.0, 0.0]
        entry.t0,           // Start time: 0
        entry.tEnd,         // End time: 2π
        0.01,               // Save interval
        1e-8                // Tolerance
    );
    
    // Compare with analytical solution
    Vector<Real> computed = solution.getSolutionAtTime(entry.tEnd);
    Vector<Real> expected = entry.endSolution;
    Real error = (computed - expected).NormL2();
    
    std::cout << "SHO error after 1 period: " << error << std::endl;
}
```

### Example 2: Test Fixed-Step Solver with Linear System

```cpp
#include "test_data/diff_eq_systems_test_bed.h"
#include "algorithms/ODESystemSolver.h"

void testFixedStepSolver() {
    // Get linear systems (simplest tests)
    auto linearSystems = ODESystemTestBed::getLinearSystems();
    
    for (const auto* entry : linearSystems) {
        // Use RK4 fixed-step solver
        ODESystemFixedStepSolver solver(
            *entry->system, 
            ODESystemStepCalculators::RK4_Basic
        );
        
        int numSteps = 1000;
        auto solution = solver.integrate(
            entry->ic, entry->t0, entry->tEnd, numSteps
        );
        
        Vector<Real> computed = solution.getSolutionAtEndpoint();
        Real error = (computed - entry->endSolution).NormL2();
        
        std::cout << entry->name << " error: " << error << std::endl;
    }
}
```

### Example 3: Test Stiff Solver with Robertson Problem

```cpp
#include "test_data/stiff_ode_test_bed.h"
#include "algorithms/ODEStiffSolvers.h"

void testStiffSolver() {
    // Get Robertson problem (extreme stiffness)
    TestStiffODE test = getRobertsonTest();
    
    // Create system with Jacobian
    auto system = test.createStiffSystem();
    
    // Solve with BDF2
    Real h = 1e-4;  // Small initial step
    auto solution = SolveBDF2(
        *system,
        test.tStart,
        test.initialCondition,
        1.0,  // Integrate to t=1 (not full 10^11)
        h,
        20,   // Max Newton iterations
        1e-10 // Newton tolerance
    );
    
    // Check mass conservation: y1 + y2 + y3 = 1
    Vector<Real> yEnd = solution.getSolutionAtEndpoint();
    Real conservation = checkConservation(yEnd, test);
    std::cout << "Mass conservation error: " << conservation << std::endl;
    
    // Check non-negativity (chemical concentrations)
    bool valid = verifyNonNegative(yEnd);
    std::cout << "Non-negative: " << (valid ? "yes" : "NO!") << std::endl;
}
```

### Example 4: Filter Systems by Category and Dimension

```cpp
void filterSystems() {
    // Get all 2D oscillatory systems
    auto oscillatory = ODESystemTestBed::getOscillatorySystems();
    auto dim2 = ODESystemTestBed::get2DSystems();
    
    std::cout << "Oscillatory systems: " << oscillatory.size() << std::endl;
    std::cout << "2D systems: " << dim2.size() << std::endl;
    
    // Get chaotic systems for long-time integration
    auto chaotic = ODESystemTestBed::getChaoticSystems();
    for (const auto* sys : chaotic) {
        std::cout << sys->name << " (dim=" << sys->system->getDim() << ")" << std::endl;
    }
}
```

### Example 5: Verify Lorenz Attractor (Chaotic System)

```cpp
#include "test_data/diff_eq_systems_test_bed.h"

void testLorenzSystem() {
    // Get Lorenz system
    auto* lorenz = ODESystemTestBed::getGeneralSystemByName("Lorenz");
    Vector<Real> ic{1.0, 1.0, 1.0};
    
    // Integrate for long time
    DormandPrince5Integrator integrator(*lorenz);
    auto solution = integrator.integrate(ic, 0.0, 100.0, 0.01, 1e-10);
    
    // Lorenz attractor should stay bounded
    Real maxNorm = 0.0;
    for (int i = 0; i < solution.getNumSteps(); ++i) {
        Real norm = solution.getSolutionAtStep(i).NormL2();
        maxNorm = std::max(maxNorm, norm);
    }
    
    // Solution should remain bounded (attractor is bounded)
    std::cout << "Max norm: " << maxNorm << " (should be < 100)" << std::endl;
}
```

### Example 6: Iterate Over All Stiff Tests by Difficulty

```cpp
void testAllStiffByDifficulty() {
    auto allTests = getAllStiffODETests();
    
    // Sort by difficulty
    std::sort(allTests.begin(), allTests.end(),
        [](const TestStiffODE& a, const TestStiffODE& b) {
            return a.difficulty < b.difficulty;
        });
    
    for (const auto& test : allTests) {
        std::cout << "Testing: " << test.name 
                  << " (dim=" << test.dimension
                  << ", stiffness=" << test.stiffnessCategory
                  << ", difficulty=" << test.difficulty << ")" << std::endl;
        
        auto system = test.createStiffSystem();
        
        // Test with appropriate solver based on stiffness
        // ... implementation depends on available solvers
    }
}
```

### Example 7: Conservation Law Verification (Kepler Problem)

```cpp
void testKeplerConservation() {
    // Get Kepler problem
    auto keplerEntry = ODESystemTestBed::getGeneralSystemByCategory(ODECategory::Physical);
    
    // Find Kepler system
    for (const auto* entry : keplerEntry) {
        if (entry->name == "Kepler") {
            // Use symplectic integrator for energy conservation
            IODESystem* sys = const_cast<IODESystem*>(entry->system);
            ODESystemLeapfrogSolver solver(*sys);
            
            Vector<Real> ic = entry->typicalIC;  // [1.0, 0.0, 0.0, 1.0]
            auto solution = solver.integrate(ic, 0.0, 100.0, 10000);
            
            // Check energy conservation
            KeplerProblemODE kepler(1.0);
            Real E0 = kepler.getTotalEnergy(ic);
            Real L0 = kepler.getAngularMomentum(ic);
            
            Vector<Real> yEnd = solution.getSolutionAtEndpoint();
            Real Ef = kepler.getTotalEnergy(yEnd);
            Real Lf = kepler.getAngularMomentum(yEnd);
            
            std::cout << "Energy drift: " << std::abs(Ef - E0) << std::endl;
            std::cout << "Angular momentum drift: " << std::abs(Lf - L0) << std::endl;
        }
    }
}
```

---

## MML ODE Solver Reference

### Explicit Methods (Non-Stiff)

| Solver | Order | Type | Use Case |
|--------|-------|------|----------|
| `ODESystemFixedStepSolver` + `Euler_StepCalculator` | 1 | Fixed | Educational |
| `ODESystemFixedStepSolver` + `RungeKutta4_StepCalculator` | 4 | Fixed | Simple problems |
| `CashKarpIntegrator` | 5(4) | Adaptive | General purpose |
| `DormandPrince5Integrator` | 5(4) | Adaptive | **Recommended default** |
| `DormandPrince8Integrator` | 8(7) | Adaptive | High precision |
| `ODESystemLeapfrogSolver` | 2 | Symplectic | Hamiltonian systems |

### Implicit Methods (Stiff)

| Solver | Order | A-Stable | Use Case |
|--------|-------|----------|----------|
| `SolveBackwardEuler()` | 1 | Yes | Very stiff, simple |
| `SolveBDF2()` | 2 | Yes | Industry standard stiff |
| `SolveRosenbrock23()` | 2/3 | Yes | Adaptive stiff |

### Solver Selection Guide

```
Is the system stiff?
├── No → Use explicit method
│   ├── Need high precision? → DormandPrince8Integrator
│   ├── Hamiltonian system? → ODESystemLeapfrogSolver
│   └── General use → DormandPrince5Integrator
│
└── Yes → Use implicit method
    ├── Extreme stiffness (ratio > 10⁸)? → SolveBackwardEuler
    ├── Moderate stiffness? → SolveBDF2
    └── Want adaptive stepping? → SolveRosenbrock23
```

---

## Test Design Patterns

### Pattern 1: Accuracy Verification
Use systems with **analytical solutions** to verify solver order:

```cpp
// Error should decrease as h^p for order-p method
for (int steps : {100, 200, 400, 800}) {
    auto sol = solver.integrate(ic, t0, tEnd, steps);
    Real error = (sol.getSolutionAtEndpoint() - exact).NormL2();
    // error_{2n} / error_n ≈ 2^{-p}
}
```

### Pattern 2: Stability Region Testing
Use **stiff systems** to test A-stability:

```cpp
// Explicit methods will fail/explode
// Implicit methods should converge
auto linearStiff = getLinearStiffTest();
// Test with h > 2/|λ_max| where explicit methods are unstable
```

### Pattern 3: Long-Time Integration
Use **conservation laws** (Kepler, Lotka-Volterra) to detect drift:

```cpp
// Symplectic methods should preserve energy
// Others will show linear drift
Real energyDrift = |E(tEnd) - E(0)|;
```

### Pattern 4: Chaotic Sensitivity
Use **chaotic systems** to test step-size adaptivity:

```cpp
// Lorenz system: small perturbations grow exponentially
// Good adaptive solver should maintain accuracy by reducing h
```

---

## See Also

- [Integration_testbed.md](Integration_testbed.md) - Numerical integration test cases
- [LinAlgSystems_testbed.md](LinAlgSystems_testbed.md) - Linear algebra test matrices
- [ODESystemSolver.h](../../mml/algorithms/ODESystemSolver.h) - Solver implementations
- [ODEAdaptiveIntegrator.h](../../mml/algorithms/ODEAdaptiveIntegrator.h) - Adaptive integrators
- [ODEStiffSolvers.h](../../mml/algorithms/ODEStiffSolvers.h) - Stiff system solvers
