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
#include "mml/base/VectorTypes.h"
#include "mml/core/CoordTransf/CoordTransfLorentz.h"
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
        0.0, 0.99, 200, "lorentz_gamma_factor.txt");

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
    
    std::ofstream ofs("results/twin_paradox_worldlines.txt");
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
        0.0, 2 * T_trip, 200, "twin_paradox_spacetime.txt");
    
    // Visualize: Proper time (aging) comparison
    LinearInterpRealFunc age_A(t_data, tau_A_data);
    LinearInterpRealFunc age_B(t_data, tau_B_data);
    
    std::vector<LinearInterpRealFunc> aging = {age_A, age_B};
    std::vector<std::string> aging_labels = {"Twin A (Earth) - biological age", "Twin B (Ship) - biological age"};
    
    Visualizer::VisualizeMultiRealFunction(aging,
        "Twin Paradox - Biological Aging Comparison",
        aging_labels,
        0.0, 2 * T_trip, 200, "twin_paradox_aging.txt");

    std::cout << "✓ Worldlines spacetime diagram complete!\n";
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
    std::cout << "  - results/twin_paradox_worldlines.txt\n";
    
    return 0;
}