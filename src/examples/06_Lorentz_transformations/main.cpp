/******************************************************************************
 * MML Example: Lorentz Transformations - Special Relativity
 * ============================================================================
 * 
 * Demonstrates MML's special relativity capabilities:
 * 
 *   1. TIME DILATION     - Moving clocks tick slower (γ factor)
 *   2. LENGTH CONTRACTION - Moving objects are shorter
 *   3. TWIN PARADOX      - The traveling twin ages less!
 * 
 * Physics: Einstein's Special Relativity (1905)
 *   - The speed of light c is constant in all inertial frames
 *   - This leads to Lorentz transformations between frames
 *   - γ = 1/√(1 - v²/c²) is the Lorentz factor
 * 
 * Units: We use c = 1 (natural units), so velocities are fractions of c.
 * 
 * Output: Numerical demonstrations of relativistic effects
 * 
 * Build: cmake --build build --target Example05_LorentzTransform
 * Run:   ./build/src/examples/Debug/Example05_LorentzTransform
 * 
 *****************************************************************************/

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "mml/base/Vector/VectorTypes.h"
#include "mml/core/CoordTransf/CoordTransfLorentz.h"
#include "mml/core/Derivation.h"
#include "mml/core/MetricTensor.h"
#include "mml/core/Integration/Integration1D.h"
#include "mml/tools/Visualizer.h"
#include "mml/tools/Serializer.h"
#endif

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>

using namespace MML;

// Calculate Lorentz factor γ
Real gamma(Real v) {
    return 1.0 / std::sqrt(1.0 - v*v);
}

/******************************************************************************
 * Helper: Proper time integrand via Minkowski metric tensor
 * 
 * For a worldline x^μ(λ) in Minkowski spacetime, the proper time is:
 *    τ = ∫ √(-g_μν dx^μ/dλ dx^ν/dλ) dλ
 * 
 * With g = diag(-1,1,1,1) and parameterizing by coordinate time t:
 *    dτ/dt = √(1 - v²) = 1/γ
 * 
 * This class computes the integrand numerically using the metric tensor
 * and numerical differentiation — no analytical formulas assumed!
 * 
 * Adapted from Chapter 17 (Minkowski spacetime path integrals).
 *****************************************************************************/
class ProperTimeIntegrand : public IRealFunction
{
    const IParametricCurve<4>& _worldline;
public:
    ProperTimeIntegrand(const IParametricCurve<4>& worldline) : _worldline(worldline) {}

    Real operator()(Real t) const override
    {
        // Numerically differentiate the 4D worldline to get tangent vector
        auto tangent = Derivation::DeriveCurve<4>(_worldline, t, nullptr);

        // Evaluate Minkowski metric: g_μν u^μ u^ν = -(dt/dλ)² + (dx/dλ)² + (dy/dλ)² + (dz/dλ)²
        MetricTensorMinkowski metric;
        Tensor2<4> g = metric(VectorN<Real, 4>{0.0, 0.0, 0.0, 0.0});

        Real ds2 = -g(tangent, tangent);   // = 1 - v²  (for timelike worldline)

        if (ds2 < 0.0)
            throw std::runtime_error("ProperTimeIntegrand: spacelike interval (v > c)!");

        return std::sqrt(ds2);
    }
};

/******************************************************************************
 * SCENARIO 1: Time Dilation - Moving Clocks Run Slow
 * 
 * A spaceship travels at various speeds. Compare elapsed time.
 *****************************************************************************/
void Demo_TimeDilation()
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  SCENARIO 1: Time Dilation - Moving Clocks Run Slow\n";
    std::cout << std::string(70, '=') << "\n\n";

    std::cout << "Imagine a spaceship traveling past Earth.\n";
    std::cout << "Earth observers measure time T (coordinate time).\n";
    std::cout << "The spaceship measures proper time τ (time on moving clock).\n\n";

    std::cout << "Relationship: τ = T / γ  where γ = 1/√(1 - v²/c²)\n\n";

    std::cout << std::fixed << std::setprecision(4);
    
    std::cout << "  Speed (v/c)    γ factor    Earth: 1 year    Ship: τ (years)\n";
    std::cout << "  -----------------------------------------------------------\n";
    
    std::vector<Real> speeds = {0.1, 0.5, 0.8, 0.9, 0.95, 0.99, 0.999};
    
    // Collect data for visualization
    std::vector<Real> v_data, gamma_data, tau_data;
    
    for (Real v : speeds) {
        Real g = gamma(v);
        Real tau = 1.0 / g;  // Proper time for 1 year Earth time
        
        v_data.push_back(v);
        gamma_data.push_back(g);
        tau_data.push_back(tau);
        
        std::cout << "     " << std::setw(5) << v 
                  << "       " << std::setw(7) << g
                  << "         1.0000          " << std::setw(7) << tau << "\n";
    }
    
    std::cout << "\n  → At 99% of c, only 0.14 years pass on the ship!\n";
    std::cout << "  → At 99.9% of c, only 0.04 years pass!\n";

    // Create visualization: Lorentz factor γ vs velocity
    // Use a lambda to create smooth curve
    auto gamma_func = [](Real v) -> Real {
        if (v >= 0.9999) return 100.0;  // Cap for display
        return 1.0 / std::sqrt(1.0 - v*v);
    };
    
    // Create interpolated function for visualization (use MML::Vector)
    int numPts = 101;
    Vector<Real> v_smooth(numPts), g_smooth(numPts);
    for (int i = 0; i < numPts; i++) {
        v_smooth[i] = i * 0.0099;  // 0 to 0.99
        g_smooth[i] = gamma_func(v_smooth[i]);
    }
    LinearInterpRealFunc gamma_curve(v_smooth, g_smooth);
    
    Visualizer::VisualizeRealFunction(gamma_curve, 
        "Lorentz Factor γ = 1/√(1-v²/c²)",
        0.0, 0.99, 200, "lorentz_gamma_factor.mml");

    std::cout << "\n✓ Time dilation demonstration complete!\n";
}


/******************************************************************************
 * SCENARIO 2: Length Contraction - Moving Objects are Shorter
 * 
 * A 100m spaceship at various speeds - how long does Earth measure it?
 *****************************************************************************/
void Demo_LengthContraction()
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  SCENARIO 2: Length Contraction\n";
    std::cout << std::string(70, '=') << "\n\n";

    std::cout << "A spaceship has proper length L₀ = 100 meters (rest length).\n";
    std::cout << "When moving, Earth observers measure contracted length L.\n\n";

    std::cout << "Relationship: L = L₀ / γ = L₀ × √(1 - v²/c²)\n\n";

    Real L0 = 100.0;  // Proper length in meters
    
    std::cout << std::fixed << std::setprecision(2);
    
    std::cout << "  Speed (v/c)    γ factor    Measured Length (m)\n";
    std::cout << "  ------------------------------------------------\n";
    
    std::vector<Real> speeds = {0.0, 0.1, 0.5, 0.8, 0.9, 0.95, 0.99};
    
    for (Real v : speeds) {
        Real g = (v == 0) ? 1.0 : gamma(v);
        Real L = L0 / g;
        
        std::cout << "     " << std::setw(4) << v 
                  << "       " << std::setw(6) << g
                  << "          " << std::setw(6) << L << "\n";
    }
    
    std::cout << "\n  → At 90% of c, the 100m ship appears only 43.6m long!\n";
    std::cout << "  → At 99% of c, it's just 14.1m!\n";
    std::cout << "\n✓ Length contraction demonstration complete!\n";
}


/******************************************************************************
 * SCENARIO 3: The Twin Paradox
 * 
 * The classic paradox: One twin stays on Earth, the other travels to a 
 * distant star and back. Who ages more?
 *****************************************************************************/
void Demo_TwinParadox()
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  SCENARIO 3: The Twin Paradox\n";
    std::cout << std::string(70, '=') << "\n\n";

    std::cout << "Setup:\n";
    std::cout << "  - Twin A stays on Earth\n";
    std::cout << "  - Twin B travels to Alpha Centauri (4 light-years) and back\n";
    std::cout << "  - Ship speed: 0.8c (80% of light speed)\n\n";

    Real distance = 4.0;   // light-years to star
    Real v = 0.8;          // speed as fraction of c
    Real g = gamma(v);
    
    // Earth perspective (coordinate time)
    Real T_outbound = distance / v;  // Time to reach star
    Real T_return = distance / v;    // Time to return
    Real T_total = T_outbound + T_return;
    
    // Ship perspective (proper time)
    Real tau_outbound = T_outbound / g;
    Real tau_return = T_return / g;
    Real tau_total = tau_outbound + tau_return;
    
    std::cout << std::fixed << std::setprecision(2);
    
    std::cout << "Journey Analysis:\n";
    std::cout << "  Distance to star: " << distance << " light-years\n";
    std::cout << "  Ship velocity:    " << v << "c\n";
    std::cout << "  Lorentz factor γ: " << g << "\n\n";
    
    std::cout << "Outbound Trip (Earth → Alpha Centauri):\n";
    std::cout << "  Earth time:  " << T_outbound << " years\n";
    std::cout << "  Ship time:   " << tau_outbound << " years\n\n";
    
    std::cout << "Return Trip (Alpha Centauri → Earth):\n";
    std::cout << "  Earth time:  " << T_return << " years\n";
    std::cout << "  Ship time:   " << tau_return << " years\n\n";
    
    std::cout << "═══════════════════════════════════════════════════════════════\n";
    std::cout << "  TOTAL JOURNEY:\n";
    std::cout << "    Twin A (Earth):  aged " << T_total << " years\n";
    std::cout << "    Twin B (Ship):   aged " << tau_total << " years\n";
    std::cout << "    Difference:      " << (T_total - tau_total) << " years younger!\n";
    std::cout << "═══════════════════════════════════════════════════════════════\n\n";
    
    std::cout << "Why is this not a paradox?\n";
    std::cout << "  The situation is NOT symmetric! Twin B accelerates and decelerates,\n";
    std::cout << "  breaking the equivalence of inertial frames. Twin B experiences\n";
    std::cout << "  non-inertial motion, so they objectively age less.\n\n";
    
    // Demonstrate with MML's Lorentz transformation
    std::cout << "Using MML's CoordTransfLorentzXAxis:\n";
    CoordTransfLorentzXAxis lorentz(v);
    
    // Event: Twin A celebrates 10th anniversary on Earth (t=10, x=0)
    Vector4Minkowski earthEvent{10.0, 0.0, 0.0, 0.0};
    Vector4Minkowski shipEvent = lorentz.transf(earthEvent);
    
    std::cout << "  Event in Earth frame: t=" << earthEvent.T() << " years, x=" << earthEvent.X() << "\n";
    std::cout << "  Same event in ship frame: t'=" << shipEvent.T() << " years, x'=" << shipEvent.X() << "\n";
    
    std::cout << "\n✓ Twin Paradox demonstration complete!\n";
}


/******************************************************************************
 * SCENARIO 4: Worldlines in Spacetime
 * 
 * Visualize the paths through spacetime (Minkowski diagram)
 *****************************************************************************/
void Demo_Worldlines()
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  SCENARIO 4: Worldlines in Spacetime\n";
    std::cout << std::string(70, '=') << "\n\n";

    std::cout << "Creating spacetime diagram (Minkowski diagram)...\n\n";
    
    Real v = 0.8;  // Ship velocity
    Real distance = 4.0;  // Distance to star (light-years)
    Real T_trip = distance / v;  // Time for one-way trip
    
    Real g = gamma(v);
    int numPoints = 100;
    Real dt = 2 * T_trip / numPoints;
    
    // Storage for visualization (use MML::Vector)
    Vector<Real> t_data(numPoints + 1);
    Vector<Real> x_A_data(numPoints + 1);  // Twin A (stays on Earth)
    Vector<Real> x_B_data(numPoints + 1);  // Twin B (traveling twin)
    Vector<Real> tau_A_data(numPoints + 1);  // Proper time A
    Vector<Real> tau_B_data(numPoints + 1);  // Proper time B
    
    Real tau_A = 0, tau_B = 0;
    
    std::ofstream ofs("results/twin_paradox_worldlines.mml");
    ofs << "# Twin Paradox Worldlines\n";
    ofs << "# t_earth  x_A(earth)  x_B(earth)  tau_A  tau_B\n";
    ofs << std::fixed << std::setprecision(4);
    
    for (int i = 0; i <= numPoints; i++) {
        Real t = i * dt;
        Real x_A = 0;  // Twin A stays at origin
        Real x_B;
        
        if (t <= T_trip) {
            x_B = v * t;  // Outbound
        } else {
            x_B = 2 * distance - v * t;  // Return
        }
        
        // Proper time accumulation
        tau_A = t;  // Earth twin: proper time = coordinate time
        if (t <= T_trip) {
            tau_B = t / g;
        } else {
            tau_B = T_trip / g + (t - T_trip) / g;
        }
        
        t_data[i] = t;
        x_A_data[i] = x_A;
        x_B_data[i] = x_B;
        tau_A_data[i] = tau_A;
        tau_B_data[i] = tau_B;
        
        ofs << t << " " << x_A << " " << x_B << " " << tau_A << " " << tau_B << "\n";
    }
    ofs.close();
    
    std::cout << "Worldline data saved to: results/twin_paradox_worldlines.txt\n";
    std::cout << "\n  Columns:\n";
    std::cout << "    t_earth     - Coordinate time (Earth frame)\n";
    std::cout << "    x_A, x_B    - Positions of Twin A and B\n";
    std::cout << "    tau_A, tau_B - Proper times (biological ages)\n\n";

    // Visualize: Spacetime diagram (worldlines)
    LinearInterpRealFunc worldline_A(t_data, x_A_data);
    LinearInterpRealFunc worldline_B(t_data, x_B_data);
    
    // Create multi-function visualization for worldlines
    std::vector<LinearInterpRealFunc> worldlines = {worldline_A, worldline_B};
    std::vector<std::string> worldline_labels = {"Twin A (stays on Earth)", "Twin B (travels to star)"};
    
    Visualizer::VisualizeMultiRealFunction(worldlines,
        "Twin Paradox - Spacetime Worldlines (position vs time)",
        worldline_labels,
        0.0, 2 * T_trip, 200, "twin_paradox_spacetime.mml");
    
    // Visualize: Proper time (aging) comparison
    LinearInterpRealFunc age_A(t_data, tau_A_data);
    LinearInterpRealFunc age_B(t_data, tau_B_data);
    
    std::vector<LinearInterpRealFunc> aging = {age_A, age_B};
    std::vector<std::string> aging_labels = {"Twin A (Earth) - biological age", "Twin B (Ship) - biological age"};
    
    Visualizer::VisualizeMultiRealFunction(aging,
        "Twin Paradox - Biological Aging Comparison",
        aging_labels,
        0.0, 2 * T_trip, 200, "twin_paradox_aging.mml");

    std::cout << "✓ Worldlines spacetime diagram complete!\n";
}


/******************************************************************************
 * SCENARIO 5: Realistic Twin Paradox — Journey to Star B
 * 
 * A proper relativistic round-trip with:
 *   Phase 1: Constant proper acceleration (1g) from rest to 0.8c
 *   Phase 2: Coast at 0.8c toward orbit entry point
 *   Phase 3: Decelerate (1g) to rest at orbit entry
 *   Phase 4: Semi-circular orbit (CW) around the star
 *   Phase 5-7: Mirror of 1-3 for return trip
 * 
 * Geometry: The ship aims for the TOP of the orbit circle (D_star, +r_orb),
 *           flies a CW semi-circle around Star B through the far side,
 *           exits at the BOTTOM (D_star, -r_orb), and returns to Earth.
 *           This creates a slight approach angle θ = atan(r_orb/D_star).
 * 
 * Physics: Hyperbolic motion (constant proper acceleration), Lorentz factor,
 *          proper time integration — all computed analytically.
 * 
 * Famous coincidence: 1g ≈ 1.032 ly/yr² (almost exactly 1 in natural units!)
 *****************************************************************************/
void Demo_Worldlines3D()
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  SCENARIO 5: Realistic Twin Paradox — Journey to Star B\n";
    std::cout << std::string(70, '=') << "\n\n";

    // ═══════════════════════════════════════════════════════════════════════
    // Physical parameters (natural units: c = 1, distances in ly, time in yr)
    // ═══════════════════════════════════════════════════════════════════════
    const Real a       = 1.032;   // 1g proper acceleration ≈ 1.032 ly/yr²
    const Real v_max   = 0.8;     // Maximum cruise velocity (0.8c)
    const Real D_star  = 10.0;    // Distance to Star B (light-years)
    const Real r_orb   = 0.5;     // Orbital radius (exaggerated for visualization)
    const Real v_orb   = 0.5;     // Orbital velocity (exaggerated for relativistic effects)
    const Real scale   = 4.0;     // Visualization scale factor

    // Derived Lorentz factors
    const Real gamma_max = gamma(v_max);   // 5/3 ≈ 1.6667
    const Real gamma_orb = gamma(v_orb);   // 1/√0.75 ≈ 1.1547

    // ═══════════════════════════════════════════════════════════════════════
    // Geometry: Ship aims for orbit entry point at top of orbit circle
    //
    //   Orbit center: (D_star, 0) = (10, 0)
    //   Entry point:  (D_star, +r_orb) — top of circle
    //   Exit point:   (D_star, -r_orb) — bottom of circle
    //   Semi-circle goes CW: top → far side → bottom
    //
    //   Approach angle: θ = atan(r_orb / D_star) ≈ 2.86°
    //   Travel distance: L = √(D_star² + r_orb²) ≈ 10.012 ly
    // ═══════════════════════════════════════════════════════════════════════
    const Real D_linear = std::sqrt(D_star*D_star + r_orb*r_orb);
    const Real cos_th   = D_star / D_linear;   // ≈ 0.9988
    const Real sin_th   = r_orb  / D_linear;   // ≈ 0.0500
    const Real theta    = std::atan2(r_orb, D_star);

    std::cout << "Mission Profile:\n";
    std::cout << "  Destination:      Star B at " << D_star << " light-years\n";
    std::cout << "  Propulsion:       Constant proper acceleration a = 1g ≈ 1.032 ly/yr²\n";
    std::cout << "  Max cruise speed: " << v_max << "c   (γ = " << std::fixed 
              << std::setprecision(4) << gamma_max << ")\n";
    std::cout << "  Orbital radius:   " << r_orb << " ly  (exaggerated for visualization)\n";
    std::cout << "  Orbital velocity: " << v_orb << "c   (γ = " << gamma_orb << ")\n";
    std::cout << "  Approach angle:   " << std::setprecision(2) << theta * 180.0 / Constants::PI 
              << "°  (aiming for orbit entry)\n";
    std::cout << "  Travel distance:  " << std::setprecision(3) << D_linear << " ly (each way)\n\n";

    // ═══════════════════════════════════════════════════════════════════════
    // Phase calculations — all analytical using hyperbolic motion equations
    //
    // Hyperbolic motion from rest with proper acceleration a:
    //   v(t) = at / √(1 + (at)²)
    //   x(t) = (√(1 + (at)²) - 1) / a    (distance along travel direction)
    //   τ(t) = arcsinh(at) / a
    // ═══════════════════════════════════════════════════════════════════════

    // ── Phase 1: Accelerate from rest to v_max along approach direction ──
    const Real at1       = v_max * gamma_max;     // = 4/3
    const Real t_accel   = at1 / a;               // Coordinate time
    const Real d_accel   = (std::sqrt(1.0 + at1*at1) - 1.0) / a;  // Distance along path
    const Real tau_accel = std::asinh(at1) / a;   // Proper time

    // ── Phase 2: Coast at v_max along approach direction ──
    const Real d_coast   = D_linear - 2.0 * d_accel;   // Remaining distance along path
    const Real t_coast   = d_coast / v_max;
    const Real tau_coast = t_coast / gamma_max;

    // ── Phase 4: Semi-circular orbit (CW from top to bottom of circle) ──
    const Real T_semi_orbit   = Constants::PI * r_orb / v_orb;     // Half the full period
    const Real tau_semi_orbit = T_semi_orbit / gamma_orb;

    // ── Phase transition times (coordinate time, Earth frame) ────────────
    const Real T1 = t_accel;                          // End of acceleration
    const Real T2 = T1 + t_coast;                     // End of coast
    const Real T3 = T2 + t_accel;                     // Arrive at orbit entry (D_star, +r_orb)
    const Real T4 = T3 + T_semi_orbit;                // Semi-orbit complete (D_star, -r_orb)
    const Real T5 = T4 + t_accel;                     // Return acceleration done
    const Real T6 = T5 + t_coast;                     // Return coast done
    const Real T7 = T6 + t_accel;                     // HOME!

    // ── Cumulative proper times at phase boundaries ──────────────────────
    const Real tau_T1 = tau_accel;
    const Real tau_T2 = tau_T1 + tau_coast;
    const Real tau_T3 = tau_T2 + tau_accel;
    const Real tau_T4 = tau_T3 + tau_semi_orbit;
    const Real tau_T5 = tau_T4 + tau_accel;
    const Real tau_T6 = tau_T5 + tau_coast;
    const Real tau_T7 = tau_T6 + tau_accel;

    // ═══════════════════════════════════════════════════════════════════════
    // Journey Analysis — Print phase-by-phase breakdown
    // ═══════════════════════════════════════════════════════════════════════
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Phase                     Earth Time   Ship Time   Distance\n";
    std::cout << "  ─────────────────────────────────────────────────────────────\n";
    std::cout << "  1. Accelerate (0→0.8c)     " << std::setw(6) << t_accel  
              << " yr    " << std::setw(6) << tau_accel << " yr   " 
              << std::setw(6) << d_accel << " ly\n";
    std::cout << "  2. Coast at 0.8c          " << std::setw(6) << t_coast  
              << " yr    " << std::setw(6) << tau_coast << " yr   " 
              << std::setw(6) << d_coast << " ly\n";
    std::cout << "  3. Decelerate (0.8c→0)     " << std::setw(6) << t_accel  
              << " yr    " << std::setw(6) << tau_accel << " yr   " 
              << std::setw(6) << d_accel << " ly\n";
    std::cout << "  4. Semi-orbit Star B      " << std::setw(6) << T_semi_orbit  
              << " yr    " << std::setw(6) << tau_semi_orbit << " yr    orbit\n";
    std::cout << "  5. Accelerate (0→0.8c)     " << std::setw(6) << t_accel  
              << " yr    " << std::setw(6) << tau_accel << " yr   "
              << std::setw(6) << d_accel << " ly\n";
    std::cout << "  6. Coast at 0.8c          " << std::setw(6) << t_coast  
              << " yr    " << std::setw(6) << tau_coast << " yr   " 
              << std::setw(6) << d_coast << " ly\n";
    std::cout << "  7. Decelerate (0.8c→0)     " << std::setw(6) << t_accel  
              << " yr    " << std::setw(6) << tau_accel << " yr   " 
              << std::setw(6) << d_accel << " ly\n";
    std::cout << "  ─────────────────────────────────────────────────────────────\n";

    std::cout << "\n";
    std::cout << "  ═══════════════════════════════════════════════════════════\n";
    std::cout << "    TOTAL JOURNEY:\n";
    std::cout << "      Twin A (Earth):   aged " << T7 << " years\n";
    std::cout << "      Twin B (Ship):    aged " << tau_T7 << " years\n";
    std::cout << "      Twin B is " << (T7 - tau_T7) << " years YOUNGER!\n";
    std::cout << "  ═══════════════════════════════════════════════════════════\n\n";

    // ═══════════════════════════════════════════════════════════════════════
    // Worldline: Twin A — stays at Earth, vertical line in spacetime
    // ═══════════════════════════════════════════════════════════════════════
    ParametricCurveFromStdFunc<3> twinA([scale](Real t) {
        return VectorN<Real, 3>{0.0, 0.0, t} * scale;
    });

    // ═══════════════════════════════════════════════════════════════════════
    // Worldline: Twin B — 7-phase realistic interstellar journey
    //
    // Outbound path aims at (D_star, +r_orb) — the top of the orbit circle.
    // d(t) is the scalar distance traveled along the approach direction.
    // Position = d(t) × (cos θ, sin θ) for outbound phases.
    //
    // CW semi-orbit: top → far side → bottom of circle
    //
    // Return path from (D_star, -r_orb) back to (0, 0).
    // Position = exit_point + d(t) × (-cos θ, +sin θ) for return phases.
    // ═══════════════════════════════════════════════════════════════════════
    ParametricCurveFromStdFunc<3> twinB(
        [=](Real t) -> VectorN<Real, 3> {
            Real x, y;

            if (t <= T1) {
                // Phase 1: Accelerate along approach direction toward (D_star, +r_orb)
                Real at = a * t;
                Real d = (std::sqrt(1.0 + at*at) - 1.0) / a;
                x = d * cos_th;
                y = d * sin_th;
            } else if (t <= T2) {
                // Phase 2: Coast at v_max along approach direction
                Real d = d_accel + v_max * (t - T1);
                x = d * cos_th;
                y = d * sin_th;
            } else if (t <= T3) {
                // Phase 3: Decelerate to rest, arriving at (D_star, +r_orb)
                Real t_rem = T3 - t;
                Real at_rem = a * t_rem;
                Real d_from_end = (std::sqrt(1.0 + at_rem*at_rem) - 1.0) / a;
                x = D_star  - d_from_end * cos_th;
                y = r_orb   - d_from_end * sin_th;
            } else if (t <= T4) {
                // Phase 4: CW semi-circle from (D_star, +r_orb) to (D_star, -r_orb)
                // angle φ goes from π/2 to -π/2 (top → right side → bottom)
                Real frac = (t - T3) / T_semi_orbit;
                Real phi = Constants::PI / 2.0 - frac * Constants::PI;
                x = D_star + r_orb * std::cos(phi);
                y = r_orb * std::sin(phi);
            } else if (t <= T5) {
                // Phase 5: Accelerate from (D_star, -r_orb) toward Earth
                // Direction: (-cos θ, +sin θ)
                Real dt = t - T4;
                Real at = a * dt;
                Real d = (std::sqrt(1.0 + at*at) - 1.0) / a;
                x = D_star - d * cos_th;
                y = -r_orb + d * sin_th;
            } else if (t <= T6) {
                // Phase 6: Coast toward Earth
                Real d = d_accel + v_max * (t - T5);
                x = D_star - d * cos_th;
                y = -r_orb + d * sin_th;
            } else {
                // Phase 7: Decelerate to rest at Earth (0, 0)
                Real t_rem = T7 - t;
                Real at_rem = a * t_rem;
                Real d_from_end = (std::sqrt(1.0 + at_rem*at_rem) - 1.0) / a;
                // Distance from origin along direction toward (D_star, -r_orb)
                x = d_from_end * cos_th;
                y = -d_from_end * sin_th;
            }

            return VectorN<Real, 3>{x, y, t} * scale;
        });

    // ── 3D Spacetime Visualization ──
    std::cout << "Creating 3D spacetime diagram...\n";
    std::cout << "  Axes: x = position (ly), y = orbital plane, z = time (yr)\n";
    std::cout << "  Outbound path angles " << std::setprecision(1) 
              << theta * 180.0 / Constants::PI << "° upward to reach orbit entry.\n";
    std::cout << "  Semi-circle at the star creates a half-helix in spacetime.\n";
    std::cout << "  Return path angles downward, arriving at Earth from below.\n\n";

    Visualizer::VisualizeMultiParamCurve3D(
        { &twinA, &twinB },
        "Realistic Twin Paradox - 3D Spacetime (x, y, t)",
        0.0, T7, 1000, "twin_paradox_realistic_3d.mml");

    // ═══════════════════════════════════════════════════════════════════════
    // Proper Time Comparison: τ_A(t) vs τ_B(t) — the aging curves
    // ═══════════════════════════════════════════════════════════════════════
    auto tauB_func = [=](Real t) -> Real {
        if (t <= T1) {
            return std::asinh(a * t) / a;
        } else if (t <= T2) {
            return tau_T1 + (t - T1) / gamma_max;
        } else if (t <= T3) {
            Real t_rem = T3 - t;
            return tau_T3 - std::asinh(a * t_rem) / a;
        } else if (t <= T4) {
            return tau_T3 + (t - T3) / gamma_orb;
        } else if (t <= T5) {
            return tau_T4 + std::asinh(a * (t - T4)) / a;
        } else if (t <= T6) {
            return tau_T5 + (t - T5) / gamma_max;
        } else {
            Real t_rem = T7 - t;
            return tau_T7 - std::asinh(a * t_rem) / a;
        }
    };

    // Sample proper time functions into vectors for visualization
    int numPts = 500;
    Vector<Real> t_pts(numPts), tau_A_pts(numPts), tau_B_pts(numPts);
    for (int i = 0; i < numPts; i++) {
        Real t = T7 * i / (numPts - 1);
        t_pts[i] = t;
        tau_A_pts[i] = t;              // Twin A: proper time = coordinate time
        tau_B_pts[i] = tauB_func(t);   // Twin B: computed analytically
    }

    LinearInterpRealFunc age_A(t_pts, tau_A_pts);
    LinearInterpRealFunc age_B(t_pts, tau_B_pts);

    std::vector<LinearInterpRealFunc> aging = {age_A, age_B};
    std::vector<std::string> aging_labels = {
        "Twin A (Earth) - biological age",
        "Twin B (Ship)  - biological age"
    };

    Visualizer::VisualizeMultiRealFunction(aging,
        "Realistic Twin Paradox - Biological Aging (proper time vs coord time)",
        aging_labels,
        0.0, T7, 500, "twin_paradox_realistic_aging.mml");

    std::cout << "  The aging comparison shows how Twin B's clock runs slower\n";
    std::cout << "  during coast phases (slope = " << std::setprecision(3) 
              << 1.0/gamma_max << " instead of 1.0)\n";
    std::cout << "  and slightly slower during orbit (slope = " 
              << 1.0/gamma_orb << ").\n\n";

    // ═══════════════════════════════════════════════════════════════════════
    // VERIFICATION: Proper Time via Minkowski Path Integral
    //
    // The "right" way to compute proper time in general relativity:
    //    τ = ∫ √(-g_μν dx^μ/dλ dx^ν/dλ) dλ
    //
    // We define Twin B's worldline as a 4D curve in Minkowski space (t,x,y,z),
    // then numerically differentiate it and apply the metric tensor.
    // This uses NO analytical proper time formulas — pure numerical integration
    // of the metric tensor along the worldline. The result should match our
    // analytical calculation to high precision.
    //
    // This demonstrates the Chapter 17 path integral approach and verifies
    // that our phase-by-phase analytical formulas are correct.
    // ═══════════════════════════════════════════════════════════════════════

    std::cout << "  ─── Verification: Minkowski Path Integral ───────────────────\n";
    std::cout << "  Computing τ = ∫ √(-g_μν dx^μ/dλ dx^ν/dλ) dλ numerically...\n\n";

    // 4D Minkowski worldline: x^μ = (t, x(t), y(t), 0)
    // Same position logic as twinB, but in 4D and WITHOUT visualization scaling
    ParametricCurveFromStdFunc<4> worldlineB(
        [=](Real t) -> VectorN<Real, 4> {
            Real x, y;

            if (t <= T1) {
                Real at = a * t;
                Real d = (std::sqrt(1.0 + at*at) - 1.0) / a;
                x = d * cos_th;
                y = d * sin_th;
            } else if (t <= T2) {
                Real d = d_accel + v_max * (t - T1);
                x = d * cos_th;
                y = d * sin_th;
            } else if (t <= T3) {
                Real t_rem = T3 - t;
                Real at_rem = a * t_rem;
                Real d_from_end = (std::sqrt(1.0 + at_rem*at_rem) - 1.0) / a;
                x = D_star  - d_from_end * cos_th;
                y = r_orb   - d_from_end * sin_th;
            } else if (t <= T4) {
                Real frac = (t - T3) / T_semi_orbit;
                Real phi = Constants::PI / 2.0 - frac * Constants::PI;
                x = D_star + r_orb * std::cos(phi);
                y = r_orb * std::sin(phi);
            } else if (t <= T5) {
                Real dt = t - T4;
                Real at = a * dt;
                Real d = (std::sqrt(1.0 + at*at) - 1.0) / a;
                x = D_star - d * cos_th;
                y = -r_orb + d * sin_th;
            } else if (t <= T6) {
                Real d = d_accel + v_max * (t - T5);
                x = D_star - d * cos_th;
                y = -r_orb + d * sin_th;
            } else {
                Real t_rem = T7 - t;
                Real at_rem = a * t_rem;
                Real d_from_end = (std::sqrt(1.0 + at_rem*at_rem) - 1.0) / a;
                x = d_from_end * cos_th;
                y = -d_from_end * sin_th;
            }

            // Minkowski 4-vector: (ct, x, y, z) with c=1
            return VectorN<Real, 4>{t, x, y, 0.0};
        });

    // Compute proper time numerically, phase by phase
    ProperTimeIntegrand integrand(worldlineB);

    // Small offsets to avoid evaluating derivatives exactly at phase boundaries
    const Real eps_t = 1e-6;

    Real tau_num_1 = IntegrateTrap(integrand, 0.0 + eps_t,      T1 - eps_t).value;
    Real tau_num_2 = IntegrateTrap(integrand, T1 + eps_t,        T2 - eps_t).value;
    Real tau_num_3 = IntegrateTrap(integrand, T2 + eps_t,        T3 - eps_t).value;
    Real tau_num_4 = IntegrateTrap(integrand, T3 + eps_t,        T4 - eps_t).value;
    Real tau_num_5 = IntegrateTrap(integrand, T4 + eps_t,        T5 - eps_t).value;
    Real tau_num_6 = IntegrateTrap(integrand, T5 + eps_t,        T6 - eps_t).value;
    Real tau_num_7 = IntegrateTrap(integrand, T6 + eps_t,        T7 - eps_t).value;
    Real tau_num_total = tau_num_1 + tau_num_2 + tau_num_3 + tau_num_4 
                       + tau_num_5 + tau_num_6 + tau_num_7;

    // Analytical proper times per phase (from our formulas)
    Real tau_ana_1 = tau_accel;
    Real tau_ana_2 = tau_coast;
    Real tau_ana_3 = tau_accel;
    Real tau_ana_4 = tau_semi_orbit;
    Real tau_ana_5 = tau_accel;
    Real tau_ana_6 = tau_coast;
    Real tau_ana_7 = tau_accel;
    Real tau_ana_total = tau_T7;

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  Phase    Analytical τ    Numerical τ     Difference\n";
    std::cout << "  ─────────────────────────────────────────────────────\n";
    auto printPhase = [](const char* name, Real ana, Real num) {
        std::cout << "  " << name << std::setw(12) << ana 
                  << " yr   " << std::setw(12) << num 
                  << " yr   " << std::setw(12) << std::abs(ana - num) << " yr\n";
    };
    printPhase("1. Accel  ", tau_ana_1, tau_num_1);
    printPhase("2. Coast  ", tau_ana_2, tau_num_2);
    printPhase("3. Decel  ", tau_ana_3, tau_num_3);
    printPhase("4. Orbit  ", tau_ana_4, tau_num_4);
    printPhase("5. Accel  ", tau_ana_5, tau_num_5);
    printPhase("6. Coast  ", tau_ana_6, tau_num_6);
    printPhase("7. Decel  ", tau_ana_7, tau_num_7);
    std::cout << "  ─────────────────────────────────────────────────────\n";
    printPhase("TOTAL     ", tau_ana_total, tau_num_total);

    Real rel_error = std::abs(tau_ana_total - tau_num_total) / tau_ana_total * 100.0;
    std::cout << "\n  Relative error: " << std::scientific << std::setprecision(2) 
              << rel_error << "%\n";
    std::cout << "  Method: τ = ∫ √(-g_μν dx^μ/dλ dx^ν/dλ) dλ  (Minkowski metric)\n";
    std::cout << "  ─────────────────────────────────────────────────────────────\n\n";

    std::cout << "✓ Realistic Twin Paradox simulation complete!\n";
}


/******************************************************************************
 * MAIN
 *****************************************************************************/
int main()
{
    std::cout << "\n";
    std::cout << "================================================================\n";
    std::cout << "     MML LORENTZ TRANSFORMATIONS - SPECIAL RELATIVITY DEMO\n";
    std::cout << "================================================================\n\n";
    
    std::cout << "   \"The distinction between past, present, and future is only\n";
    std::cout << "    a stubbornly persistent illusion.\" - Albert Einstein\n";
    std::cout << "\n";
    std::cout << "   Using natural units where c = 1 (speed of light).\n";
    std::cout << "   Time in years, distance in light-years.\n";
    std::cout << "================================================================\n";
    
    Demo_TimeDilation();
    Demo_LengthContraction();
    Demo_TwinParadox();
    Demo_Worldlines();
    Demo_Worldlines3D();
    
    std::cout << "\n";
    std::cout << "================================================================\n";
    std::cout << "         SPECIAL RELATIVITY DEMONSTRATIONS COMPLETE!\n";
    std::cout << "================================================================\n\n";
    
    std::cout << "Key Takeaways:\n";
    std::cout << "  1. Moving clocks run SLOW (time dilation)\n";
    std::cout << "  2. Moving objects are SHORTER (length contraction)\n";
    std::cout << "  3. The traveling twin ages LESS (twin paradox)\n";
    std::cout << "  4. These are REAL effects, not optical illusions!\n\n";
    
    std::cout << "Output files:\n";
    std::cout << "  - results/twin_paradox_worldlines.mml\n";
    std::cout << "  - results/twin_paradox_realistic_3d.mml\n";
    std::cout << "  - results/twin_paradox_realistic_aging.mml\n";
    
    return 0;
}