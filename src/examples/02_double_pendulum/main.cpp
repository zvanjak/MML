/******************************************************************************
 * MML Example: Double Pendulum - Chaos Demonstration
 * ============================================================================
 * 
 * Demonstrates deterministic chaos through the double pendulum:
 * 
 *   1. SINGLE TRAJECTORY - Watch the beautiful chaotic motion
 *   2. BUTTERFLY EFFECT  - Tiny initial difference → huge divergence!
 *   3. ENERGY CHECK      - Verify conservation (numerical accuracy)
 * 
 * Physics: Two coupled nonlinear differential equations
 *   - Lagrangian mechanics gives equations of motion
 *   - Highly sensitive to initial conditions (chaos!)
 *   - Total energy should be conserved
 * 
 * Output: Trajectory visualizations using MML's Visualizer
 * 
 * Build: cmake --build build --target Example02_DoublePendulum
 * Run:   ./build/src/examples/Debug/Example02_DoublePendulum
 * 
 *****************************************************************************/

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "mml/base/BaseUtils.h"
#include "mml/core/Derivation.h"
#include "mml/algorithms/ODESolvers/ODEAdaptiveIntegrator.h"
#include "mml/algorithms/ODESolvers/ODEFixedStepIntegrators.h"
#include "mml/algorithms/ODESolvers/ODESystemStepCalculators.h"
#include "mml/tools/Visualizer.h"
#include "mml/tools/Serializer.h"
#endif

// Self-contained double pendulum physics (no MPL dependency)
#include "DoublePendulum.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace MML;
using namespace Pendulum;

// Helper: Calculate total energy of double pendulum
Real calcTotalEnergy(Real m1, Real m2, Real l1, Real l2, 
                     Real theta1, Real omega1, Real theta2, Real omega2)
{
    const Real g = 9.81;
    
    // Kinetic energy
    Real T = 0.5 * m1 * l1 * l1 * omega1 * omega1
           + 0.5 * m2 * (l1 * l1 * omega1 * omega1 
                       + l2 * l2 * omega2 * omega2 
                       + 2 * l1 * l2 * omega1 * omega2 * cos(theta1 - theta2));
    
    // Potential energy (relative to pivot)
    Real V = -(m1 + m2) * g * l1 * cos(theta1) - m2 * g * l2 * cos(theta2);
    
    return T + V;
}

// Helper: Convert angles to Cartesian coordinates for visualization
void anglesToCartesian(Real l1, Real l2, Real theta1, Real theta2,
                       Real& x1, Real& y1, Real& x2, Real& y2)
{
    x1 = l1 * sin(theta1);
    y1 = -l1 * cos(theta1);
    x2 = x1 + l2 * sin(theta2);
    y2 = y1 - l2 * cos(theta2);
}

/// @brief Wrap angle to [-π, π] for visualization
static Real wrapAngle(Real theta) {
    theta = std::fmod(theta + Constants::PI, 2.0 * Constants::PI);
    if (theta < 0) theta += 2.0 * Constants::PI;
    return theta - Constants::PI;
}

/// @brief Create wrapped-angle vectors from solution data
static Vector<Real> wrapAngles(const Vector<Real>& angles) {
    Vector<Real> wrapped(angles.size());
    for (int i = 0; i < angles.size(); i++)
        wrapped[i] = wrapAngle(angles[i]);
    return wrapped;
}


/******************************************************************************
 * SCENARIO 1: Single Trajectory - Chaotic Motion
 * 
 * Start from a dramatic initial position and watch the chaos unfold!
 *****************************************************************************/
void Demo_SingleTrajectory()
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  SCENARIO 1: Chaotic Double Pendulum Motion\n";
    std::cout << std::string(70, '=') << "\n\n";

    // Physical parameters
    Real m1 = 1.0, m2 = 1.0;   // masses (kg)
    Real l1 = 1.0, l2 = 1.0;   // lengths (m)

    // Initial conditions - start from dramatic position
    Real theta1_init = Utils::DegToRad(135.0);  // First arm at 135°
    Real theta2_init = Utils::DegToRad(135.0);  // Second arm at 135°
    Real omega1_init = 0.0;                      // No initial angular velocity
    Real omega2_init = 0.0;

    std::cout << "Physical Parameters:\n";
    std::cout << "  Mass 1: " << m1 << " kg, Length 1: " << l1 << " m\n";
    std::cout << "  Mass 2: " << m2 << " kg, Length 2: " << l2 << " m\n\n";

    std::cout << "Initial Conditions:\n";
    std::cout << "  Theta 1: " << Utils::RadToDeg(theta1_init) << " degrees\n";
    std::cout << "  Theta 2: " << Utils::RadToDeg(theta2_init) << " degrees\n";
    std::cout << "  (Both arms pointing up-left - dramatic start!)\n\n";

    // Create ODE system
    DoublePendulumODE pendulum(m1, m2, l1, l2);

    // Solve with fixed-step RK4 (fast!)
    Real t_end = 30.0;  // 30 seconds (longer for richer phase space)
    int numSteps = 6000;
    Vector<Real> initCond{ theta1_init, omega1_init, theta2_init, omega2_init };

    std::cout << "Solving for " << t_end << " seconds... " << std::flush;
    
    ODESystemFixedStepSolver solver(pendulum, StepCalculators::RK4_Basic);
    ODESystemSolution sol = solver.integrate(initCond, 0.0, t_end, numSteps);
    
    std::cout << "Done! (" << sol.getTValues().size() << " time steps)\n\n";

    // Get solutions
    Vector<Real> t_vals = sol.getTValues();
    Vector<Real> theta1_vals = sol.getXValues(0);
    Vector<Real> omega1_vals = sol.getXValues(1);
    Vector<Real> theta2_vals = sol.getXValues(2);
    Vector<Real> omega2_vals = sol.getXValues(3);

    // Energy check
    Real E_init = calcTotalEnergy(m1, m2, l1, l2, theta1_init, omega1_init, theta2_init, omega2_init);
    Real E_final = calcTotalEnergy(m1, m2, l1, l2, 
                                   theta1_vals[theta1_vals.size()-1], omega1_vals[omega1_vals.size()-1],
                                   theta2_vals[theta2_vals.size()-1], omega2_vals[omega2_vals.size()-1]);
    
    std::cout << "Energy Conservation Check:\n";
    std::cout << "  Initial energy: " << std::fixed << std::setprecision(6) << E_init << " J\n";
    std::cout << "  Final energy:   " << E_final << " J\n";
    std::cout << "  Energy drift:   " << std::abs((E_final - E_init) / E_init) * 100 << "%\n\n";

    // Visualize angles over time (wrapped to [-π, π])
    Vector<Real> theta1_wrapped = wrapAngles(theta1_vals);
    Vector<Real> theta2_wrapped = wrapAngles(theta2_vals);
    PolynomInterpRealFunc theta1_func(t_vals, theta1_wrapped, 3);
    PolynomInterpRealFunc theta2_func(t_vals, theta2_wrapped, 3);

    Visualizer::VisualizeMultiRealFunction(
        std::vector<IRealFunction*>{ &theta1_func, &theta2_func },
        "Double Pendulum Angles vs Time (wrapped)",
        { "Theta 1", "Theta 2" },
        0.0, t_end, 1000, "double_pendulum_angles.mml");

    // Phase space trajectory (theta1 vs theta2) - already wrapped

    Matrix<Real> phase_points(t_vals.size(), 2);
    for (int i = 0; i < t_vals.size(); i++)
    {
        phase_points(i, 0) = wrapAngle(theta1_vals[i]);
        phase_points(i, 1) = wrapAngle(theta2_vals[i]);
    }
    SplineInterpParametricCurve<2> phase_curve(0.0, 1.0, phase_points);
    
    Visualizer::VisualizeParamCurve2D(phase_curve, "Double Pendulum Phase Space",
        0.0, 1.0, t_vals.size(), "double_pendulum_phase_space.mml");

    // Output trajectory file for animation
    std::ofstream ofs("results/double_pendulum_trajectory.mml");
    ofs << "# Double Pendulum Trajectory\n";
    ofs << "# t theta1 theta2 x1 y1 x2 y2\n";
    ofs << std::fixed << std::setprecision(6);
    
    for (int i = 0; i < t_vals.size(); i += 5)  // Every 5th point
    {
        Real x1, y1, x2, y2;
        anglesToCartesian(l1, l2, theta1_vals[i], theta2_vals[i], x1, y1, x2, y2);
        ofs << t_vals[i] << " " << theta1_vals[i] << " " << theta2_vals[i] 
            << " " << x1 << " " << y1 << " " << x2 << " " << y2 << "\n";
    }
    ofs.close();

    std::cout << "Visualizations saved to:\n";
    std::cout << "  - results/double_pendulum_angles.txt\n";
    std::cout << "  - results/double_pendulum_phase_space.txt\n";
    std::cout << "  - results/double_pendulum_trajectory.txt\n";
    std::cout << "\n✓ Single trajectory complete!\n";
}


/******************************************************************************
 * SCENARIO 2: Butterfly Effect - Chaos Demonstration
 * 
 * Two pendulums with TINY initial difference (0.001°) → completely different!
 * This is the essence of deterministic chaos.
 *****************************************************************************/
void Demo_ButterflyEffect()
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  SCENARIO 2: Butterfly Effect - Chaos Demonstration\n";
    std::cout << std::string(70, '=') << "\n\n";

    Real m1 = 1.0, m2 = 1.0;
    Real l1 = 1.0, l2 = 1.0;

    // Two nearly identical initial conditions - use high-energy start for chaos
    Real theta1_A = Utils::DegToRad(135.0);
    Real theta2_A = Utils::DegToRad(135.0);
    
    Real epsilon = Utils::DegToRad(0.001);  // Just 0.001 degree difference!
    Real theta1_B = theta1_A + epsilon;
    Real theta2_B = theta2_A;

    std::cout << "Initial Conditions:\n";
    std::cout << "  Pendulum A: theta1 = 135.000°, theta2 = 135°\n";
    std::cout << "  Pendulum B: theta1 = 135.001°, theta2 = 135°\n";
    std::cout << "  Difference: just 0.001 degrees!\n\n";

    // Create and solve both systems with fixed-step RK4
    DoublePendulumODE pendulum(m1, m2, l1, l2);
    ODESystemFixedStepSolver solver(pendulum, StepCalculators::RK4_Basic);

    Real t_end = 30.0;  // 30 seconds to see divergence
    int numSteps = 6000;
    
    Vector<Real> initCond_A{ theta1_A, 0.0, theta2_A, 0.0 };
    Vector<Real> initCond_B{ theta1_B, 0.0, theta2_B, 0.0 };

    std::cout << "Solving both pendulums... " << std::flush;
    
    ODESystemSolution sol_A = solver.integrate(initCond_A, 0.0, t_end, numSteps);
    ODESystemSolution sol_B = solver.integrate(initCond_B, 0.0, t_end, numSteps);
    
    std::cout << "Done!\n\n";

    // Get solutions
    Vector<Real> t_vals = sol_A.getTValues();
    Vector<Real> theta1_A_vals = sol_A.getXValues(0);
    Vector<Real> theta1_B_vals = sol_B.getXValues(0);

    // Calculate divergence over time
    std::cout << "Divergence over time:\n";
    std::cout << std::setw(10) << "Time (s)" << std::setw(20) << "Angle Difference\n";
    std::cout << std::string(30, '-') << "\n";
    
    std::vector<Real> divergence_times{ 0, 2, 5, 8, 12, 16, 20, 25, 30 };
    for (Real check_t : divergence_times)
    {
        // Find closest time index
        int idx = 0;
        for (int i = 0; i < t_vals.size(); i++)
        {
            if (t_vals[i] >= check_t) { idx = i; break; }
        }
        
        Real diff = std::abs(theta1_A_vals[idx] - theta1_B_vals[idx]);
        std::cout << std::fixed << std::setprecision(1) << std::setw(10) << check_t 
                  << std::setprecision(4) << std::setw(18) << Utils::RadToDeg(diff) << " deg\n";
    }
    std::cout << "\n  Initial difference: 0.001°\n";
    std::cout << "  → Grows EXPONENTIALLY over time!\n";
    std::cout << "  This is DETERMINISTIC CHAOS!\n\n";

    // Visualize both trajectories (wrapped to [-π, π])
    Vector<Real> theta1_A_wrapped = wrapAngles(theta1_A_vals);
    Vector<Real> theta1_B_wrapped = wrapAngles(theta1_B_vals);
    PolynomInterpRealFunc theta1_A_func(t_vals, theta1_A_wrapped, 3);
    PolynomInterpRealFunc theta1_B_func(t_vals, theta1_B_wrapped, 3);

    Visualizer::VisualizeMultiRealFunction(
        std::vector<IRealFunction*>{ &theta1_A_func, &theta1_B_func },
        "Butterfly Effect: Pendulum A vs B (theta1, wrapped)",
        { "Pendulum A", "Pendulum B (+ 0.001 deg)" },
        0.0, t_end, 1000, "double_pendulum_butterfly.mml");

    std::cout << "Visualization saved to: results/double_pendulum_butterfly.txt\n";
    std::cout << "\n✓ Butterfly effect demonstration complete!\n";
}


/******************************************************************************
 * MAIN
 *****************************************************************************/
int main()
{
    std::cout << "\n";
    std::cout << "======================================================================\n";
    std::cout << "         MML DOUBLE PENDULUM - CHAOS DEMONSTRATION\n";
    std::cout << "         ==========================================\n";
    std::cout << "\n";
    std::cout << "   The double pendulum is a classic example of deterministic chaos.\n";
    std::cout << "   The motion is completely determined by the equations, yet\n";
    std::cout << "   impossible to predict long-term due to extreme sensitivity\n";
    std::cout << "   to initial conditions.\n";
    std::cout << "\n";
    std::cout << "   \"Does the flap of a butterfly's wings in Brazil set off\n";
    std::cout << "    a tornado in Texas?\" - Edward Lorenz, 1972\n";
    std::cout << "======================================================================\n";

    Demo_SingleTrajectory();
    Demo_ButterflyEffect();

    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  DOUBLE PENDULUM SIMULATIONS COMPLETE!\n";
    std::cout << std::string(70, '=') << "\n\n";

    std::cout << "Output files in results/:\n";
    std::cout << "  - double_pendulum_angles.txt\n";
    std::cout << "  - double_pendulum_phase_space.txt\n";
    std::cout << "  - double_pendulum_trajectory.txt\n";
    std::cout << "  - double_pendulum_butterfly.txt\n\n";

    std::cout << "MML Visualization:\n";
    std::cout << "  All visualizations created using MML's Visualizer class.\n";
    std::cout << "  Phase space plot shows the beautiful chaotic attractor!\n\n";

    return 0;
}
