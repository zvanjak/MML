/******************************************************************************
 * MML Example: Projectile Launch with Air Resistance
 * ============================================================================
 * 
 * Demonstrates MML's projectile motion simulation capabilities:
 * 
 *   1. VACUUM VS AIR  - How air resistance drastically changes trajectories
 *   2. OPTIMAL ANGLE  - Finding the best launch angle (spoiler: not 45°!)
 *   3. BASEBALL TYPES - Different ball surfaces, different drag
 * 
 * Physics: F = mg - kv² (quadratic drag)
 *   - Drag coefficient depends on speed, shape, surface roughness
 *   - Air density decreases with altitude (isothermal/adiabatic models)
 * 
 * Output: Trajectory visualizations using MML's Visualizer
 * 
 * Build: cmake --build build --target Example01_ProjectileLaunch
 * Run:   ./build/src/examples/Debug/Example01_ProjectileLaunch
 * 
 *****************************************************************************/

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "mml/base/BaseUtils.h"
#include "mml/tools/Visualizer.h"
#include "mml/tools/Serializer.h"
#include "mml/algorithms/ODESystemSolver.h"
#endif

// Self-contained projectile physics (no MPL dependency)
#include "Projectiles2D.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace MML;
using namespace Projectile;

/******************************************************************************
 * SCENARIO 1: Vacuum vs Air Resistance Comparison
 * 
 * The classic demonstration: how much does air resistance matter?
 * Answer: A LOT! Range can be reduced by 60% or more.
 *****************************************************************************/
void Demo_VacuumVsAir()
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  SCENARIO 1: Vacuum vs Air Resistance\n";
    std::cout << std::string(70, '=') << "\n\n";

    // Launch parameters
    Real velocity = 100.0;      // m/s (like a fast baseball pitch or artillery)
    Real angle = Utils::DegToRad(45.0);  // 45 degrees
    Real initHeight = 0.0;      // ground level

    std::cout << "Launch Parameters:\n";
    std::cout << "  Initial velocity: " << velocity << " m/s\n";
    std::cout << "  Launch angle: 45 degrees\n";
    std::cout << "  Initial height: " << initHeight << " m\n\n";

    // Create ODE systems
    Projectile2DInVacuumODE vacuumSys;
    Projectile2DWithAirResistanceODE airSys(0.01);  // drag coefficient

    // Create solvers
    ProjectileMotionSolver2D vacuumSolver(vacuumSys);
    ProjectileMotionSolver2D airSolver(airSys);

    // Calculate vacuum trajectory analytically
    Real vacuumRange = vacuumSys.CalcRange(angle, initHeight, velocity);
    Real vacuumTime = vacuumSys.TimeOfFlight(angle, initHeight, velocity);
    Real vacuumMaxHeight = vacuumSys.MaxHeight(angle, initHeight, velocity);

    std::cout << "Vacuum (analytical):\n";
    std::cout << "  Range: " << std::fixed << std::setprecision(2) << vacuumRange << " m\n";
    std::cout << "  Time of flight: " << vacuumTime << " s\n";
    std::cout << "  Max height: " << vacuumMaxHeight << " m\n\n";

    // Solve numerically
    Real tMax = 1.2 * vacuumTime;  // a bit longer than vacuum time
    Real dt = 0.1;

    ProjectileTrajectory2D vacuumResult = vacuumSolver.solveEuler(angle, initHeight, velocity, tMax, dt);
    ProjectileTrajectory2D airResult = airSolver.solveEuler(angle, initHeight, velocity, tMax, dt);

    std::cout << "With Air Resistance (numerical):\n";
    std::cout << "  Range: " << airResult._range << " m\n";
    std::cout << "  Time of flight: " << airResult._timeOfFlight << " s\n";
    std::cout << "  Range reduction: " << (1 - airResult._range / vacuumRange) * 100 << "%\n\n";

    // Visualize comparison
    std::vector<LinearInterpRealFunc> trajectories;
    trajectories.push_back(vacuumResult.getYOfX());
    trajectories.push_back(airResult.getYOfX());

    std::vector<std::string> labels{"Vacuum", "With Air Resistance"};

    Visualizer::VisualizeMultiRealFunction(trajectories,
        "Projectile: Vacuum vs Air Resistance (v=" + std::to_string((int)velocity) + " m/s)",
        labels, 0, 1.1 * vacuumRange, 500, "projectile_vacuum_vs_air.txt");

    std::cout << "Visualization saved to: results/projectile_vacuum_vs_air.txt\n";
    std::cout << "\n✓ Vacuum vs Air comparison complete!\n";
}


/******************************************************************************
 * SCENARIO 2: Optimal Launch Angle
 * 
 * In vacuum, 45° gives maximum range. With air resistance?
 * The optimal angle is LOWER (typically 35-40°)!
 *****************************************************************************/
void Demo_OptimalAngle()
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  SCENARIO 2: Optimal Launch Angle\n";
    std::cout << std::string(70, '=') << "\n\n";

    Real velocity = 50.0;       // m/s
    Real initHeight = 0.0;

    // Test angles from 20° to 70°
    Vector<Real> angles(11);
    for (int i = 0; i < 11; i++)
        angles[i] = Utils::DegToRad(20.0 + i * 5.0);  // 20, 25, 30, ... 70 degrees

    std::cout << "Testing angles: 20° to 70° (step 5°)\n";
    std::cout << "Initial velocity: " << velocity << " m/s\n\n";

    // With air resistance
    Projectile2DWithAirResistanceODE airSys(0.01);
    ProjectileMotionSolver2D airSolver(airSys);

    Vector<ProjectileTrajectory2D> results = airSolver.solveForAnglesEuler(angles, initHeight, velocity, 0.1);

    // Find optimal angle
    Real maxRange = 0;
    Real optimalAngle = 0;
    
    std::cout << "Results with air resistance:\n";
    std::cout << std::setw(10) << "Angle" << std::setw(15) << "Range (m)" << "\n";
    std::cout << std::string(25, '-') << "\n";
    
    for (int i = 0; i < results.size(); i++)
    {
        Real angleDeg = Utils::RadToDeg(results[i]._angle);
        std::cout << std::fixed << std::setprecision(1);
        std::cout << std::setw(10) << angleDeg << "°" 
                  << std::setw(14) << results[i]._range << "\n";
        
        if (results[i]._range > maxRange)
        {
            maxRange = results[i]._range;
            optimalAngle = angleDeg;
        }
    }

    std::cout << "\n  OPTIMAL ANGLE: " << optimalAngle << "° (range: " << maxRange << " m)\n";
    std::cout << "  (Compare to vacuum optimal of 45°!)\n\n";

    // Visualize all trajectories
    std::vector<LinearInterpRealFunc> trajectories;
    std::vector<std::string> labels;
    
    for (auto& result : results)
    {
        trajectories.push_back(result.getYOfX());
        labels.push_back(std::to_string((int)Utils::RadToDeg(result._angle)) + " deg");
    }

    Visualizer::VisualizeMultiRealFunction(trajectories,
        "Projectile Trajectories at Different Angles",
        labels, 0, 1.1 * maxRange, 500, "projectile_angles.txt");

    std::cout << "Visualization saved to: results/projectile_angles.txt\n";
    std::cout << "\n✓ Optimal angle analysis complete!\n";
}


/******************************************************************************
 * SCENARIO 3: Baseball Comparison
 * 
 * Different ball surfaces = different drag coefficients!
 * Smooth, normal, and rough baseballs behave differently.
 *****************************************************************************/
void Demo_BaseballTypes()
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  SCENARIO 3: Baseball Surface Comparison\n";
    std::cout << std::string(70, '=') << "\n\n";

    Real velocity = 45.0;       // m/s (~100 mph pitch)
    Real angle = Utils::DegToRad(35.0);  // typical baseball hit angle
    Real initHeight = 1.0;      // bat height

    std::cout << "Baseball hit parameters:\n";
    std::cout << "  Exit velocity: " << velocity << " m/s (~100 mph)\n";
    std::cout << "  Launch angle: 35 degrees\n";
    std::cout << "  Initial height: " << initHeight << " m\n\n";

    // Different baseball types
    BaseballWithDragCoeffDependentOnSpeedODE smoothBall(BaseballType::Smooth);
    BaseballWithDragCoeffDependentOnSpeedODE normalBall(BaseballType::Normal);
    BaseballWithDragCoeffDependentOnSpeedODE roughBall(BaseballType::Rough);

    ProjectileMotionSolver2D smoothSolver(smoothBall);
    ProjectileMotionSolver2D normalSolver(normalBall);
    ProjectileMotionSolver2D roughSolver(roughBall);

    Real tMax = 10.0;
    Real dt = 0.05;

    ProjectileTrajectory2D smoothResult = smoothSolver.solveEuler(angle, initHeight, velocity, tMax, dt);
    ProjectileTrajectory2D normalResult = normalSolver.solveEuler(angle, initHeight, velocity, tMax, dt);
    ProjectileTrajectory2D roughResult = roughSolver.solveEuler(angle, initHeight, velocity, tMax, dt);

    std::cout << "Results by ball type:\n";
    std::cout << std::setw(12) << "Type" << std::setw(15) << "Range (m)" << "\n";
    std::cout << std::string(27, '-') << "\n";
    std::cout << std::setw(12) << "Smooth" << std::setw(15) << std::fixed << std::setprecision(2) << smoothResult._range << "\n";
    std::cout << std::setw(12) << "Normal" << std::setw(15) << normalResult._range << "\n";
    std::cout << std::setw(12) << "Rough" << std::setw(15) << roughResult._range << "\n\n";

    // Visualize
    std::vector<LinearInterpRealFunc> trajectories;
    trajectories.push_back(smoothResult.getYOfX());
    trajectories.push_back(normalResult.getYOfX());
    trajectories.push_back(roughResult.getYOfX());

    std::vector<std::string> labels{"Smooth ball", "Normal ball", "Rough ball"};

    Real maxRange = std::max({smoothResult._range, normalResult._range, roughResult._range});
    
    Visualizer::VisualizeMultiRealFunction(trajectories,
        "Baseball Trajectories by Surface Type",
        labels, 0, 1.1 * maxRange, 500, "projectile_baseball_types.txt");

    std::cout << "Visualization saved to: results/projectile_baseball_types.txt\n";
    std::cout << "\n✓ Baseball comparison complete!\n";
}


/******************************************************************************
 * MAIN
 *****************************************************************************/
int main()
{
    std::cout << "\n";
    std::cout << "======================================================================\n";
    std::cout << "         MML PROJECTILE LAUNCH SIMULATOR\n";
    std::cout << "         ================================\n";
    std::cout << "\n";
    std::cout << "   Demonstrating projectile motion with various drag models:\n";
    std::cout << "   - Vacuum (ideal parabola)\n";
    std::cout << "   - Constant air resistance\n";
    std::cout << "   - Speed-dependent drag (baseball physics)\n";
    std::cout << "\n";
    std::cout << "   Physics: F = mg - kv^2 (quadratic drag)\n";
    std::cout << "======================================================================\n";

    // Run all demos
    Demo_VacuumVsAir();
    Demo_OptimalAngle();
    Demo_BaseballTypes();

    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  ALL PROJECTILE SIMULATIONS COMPLETE!\n";
    std::cout << std::string(70, '=') << "\n\n";

    std::cout << "Output files created in results/:\n";
    std::cout << "  - projectile_vacuum_vs_air.txt\n";
    std::cout << "  - projectile_angles.txt\n";
    std::cout << "  - projectile_baseball_types.txt\n\n";

    std::cout << "MML Visualization:\n";
    std::cout << "  Trajectories automatically visualized using MML's Visualizer\n";
    std::cout << "  with gnuplot integration for publication-quality plots.\n\n";

    return 0;
}
