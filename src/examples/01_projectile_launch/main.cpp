/******************************************************************************
 * MML Example: Projectile Launch with Air Resistance
 * ============================================================================
 * 
 * Demonstrates MML's projectile motion simulation capabilities using
 * adaptive ODE integration with event detection for ground impact:
 * 
 *   1. VACUUM VS AIR   - Trajectory comparison + solver precision showdown
 *   2. OPTIMAL ANGLE   - Finding the best launch angle (spoiler: not 45°!)
 *   3. BASEBALL TYPES  - Different ball surfaces, different drag
 *   4. ALTITUDE DRAG   - How changing air density affects long-range shots
 * 
 * Physics: F = mg - kv² (quadratic drag)
 *   - Drag coefficient depends on speed, shape, surface roughness
 *   - Air density decreases with altitude (isothermal/adiabatic models)
 * 
 * Key feature: IODESystemWithEvents + DormandPrince5 adaptive integrator
 * detects ground impact via zero-crossing of y(t) to ~1e-12 precision.
 * No more "if(y<0) zero derivatives" hacks!
 * 
 * Build: cmake --build build --target Example01_ProjectileLaunch
 * Run:   ./build/src/examples/Release/Example01_ProjectileLaunch
 * 
 *****************************************************************************/

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "mml/base/BaseUtils.h"
#include "mml/tools/Visualizer.h"
#include "mml/tools/Serializer.h"
#include "mml/algorithms/ODESolvers/ODEAdaptiveIntegrator.h"
#endif

// Self-contained projectile physics
#include "Projectiles2D.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace MML;
using namespace Projectile;

/******************************************************************************
 * SCENARIO 1: Vacuum vs Air Resistance + Solver Precision Showdown
 * 
 * The classic demonstration: how much does air resistance matter?
 * Now with a twist: compare Euler (dt=0.1) vs DormandPrince5 with event
 * detection against the exact analytical solution.
 *****************************************************************************/
void Demo_VacuumVsAir()
{
    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  SCENARIO 1: Vacuum vs Air Resistance + Solver Precision         ║\n";
    std::cout << "╠═══════════════════════════════════════════════════════════════════╣\n";
    std::cout << "║  Comparing Euler (dt=0.1) vs DormandPrince5 + event detection    ║\n";
    std::cout << "║  against exact analytical solution for vacuum trajectory          ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n\n";

    Real velocity = 100.0;
    Real angle = Utils::DegToRad(45.0);
    Real initHeight = 0.0;

    std::cout << "Launch: v0 = " << velocity << " m/s, angle = 45 deg, h0 = " << initHeight << " m\n\n";

    // ── Analytical solution (exact) ──
    Projectile2DInVacuumODE vacuumSys;
    Real exactRange = vacuumSys.CalcRange(angle, initHeight, velocity);
    Real exactTime  = vacuumSys.TimeOfFlight(angle, initHeight, velocity);
    Real exactMaxH  = vacuumSys.MaxHeight(angle, initHeight, velocity);

    std::cout << "EXACT (analytical):\n";
    std::cout << "  Range:          " << std::fixed << std::setprecision(6) << exactRange << " m\n";
    std::cout << "  Time of flight: " << exactTime << " s\n";
    std::cout << "  Max height:     " << exactMaxH << " m\n\n";

    // ── Euler solver (dt=0.1) ──
    ProjectileMotionSolver2D vacuumSolverEuler(vacuumSys);
    Real tMax = 1.2 * exactTime;
    ProjectileTrajectory2D eulerResult = vacuumSolverEuler.solveEuler(angle, initHeight, velocity, tMax, 0.1);

    std::cout << "EULER (dt = 0.1s, " << eulerResult._tValues.size() << " points):\n";
    std::cout << "  Range:          " << eulerResult._range << " m\n";
    std::cout << "  Time of flight: " << eulerResult._timeOfFlight << " s\n";
    std::cout << "  Range error:    " << std::scientific << std::setprecision(2)
              << std::abs(eulerResult._range - exactRange) << " m ("
              << std::fixed << std::setprecision(4)
              << 100.0 * std::abs(eulerResult._range - exactRange) / exactRange << "%)\n";
    std::cout << "  Time error:     " << std::scientific << std::setprecision(2)
              << std::abs(eulerResult._timeOfFlight - exactTime) << " s\n\n";

    // ── DormandPrince5 + event detection ──
    ProjectileMotionSolver2D vacuumSolverDP5(vacuumSys);
    ProjectileTrajectory2D dp5Result = vacuumSolverDP5.solveAdaptive(angle, initHeight, velocity);

    std::cout << "DORMAND-PRINCE 5 + EVENT DETECTION (" << dp5Result._tValues.size() << " points):\n";
    std::cout << "  Range:          " << std::fixed << std::setprecision(6) << dp5Result._range << " m\n";
    std::cout << "  Time of flight: " << dp5Result._timeOfFlight << " s\n";
    std::cout << "  Range error:    " << std::scientific << std::setprecision(2)
              << std::abs(dp5Result._range - exactRange) << " m ("
              << std::fixed << std::setprecision(10)
              << 100.0 * std::abs(dp5Result._range - exactRange) / exactRange << "%)\n";
    std::cout << "  Time error:     " << std::scientific << std::setprecision(2)
              << std::abs(dp5Result._timeOfFlight - exactTime) << " s\n\n";

    // ── Precision improvement factor ──
    Real eulerRangeErr = std::abs(eulerResult._range - exactRange);
    Real dp5RangeErr   = std::abs(dp5Result._range - exactRange);
    if (dp5RangeErr > 0)
    {
        std::cout << "  >> DormandPrince5 is " << std::fixed << std::setprecision(0) 
                  << eulerRangeErr / dp5RangeErr << "x more precise than Euler!\n\n";
    }

    // ── Air resistance comparison (using adaptive solver) ──
    Projectile2DWithAirResistanceODE airSys(0.01);
    ProjectileMotionSolver2D airSolver(airSys);
    ProjectileTrajectory2D airResult = airSolver.solveAdaptive(angle, initHeight, velocity);

    std::cout << "WITH AIR RESISTANCE (k=0.01, adaptive + events):\n";
    std::cout << "  Range:          " << std::fixed << std::setprecision(2) << airResult._range << " m\n";
    std::cout << "  Time of flight: " << airResult._timeOfFlight << " s\n";
    std::cout << "  Range reduction: " << std::setprecision(1)
              << (1 - airResult._range / exactRange) * 100 << "% vs vacuum\n\n";

    // ── Visualize ──
    std::vector<LinearInterpRealFunc> trajectories;
    trajectories.push_back(dp5Result.getYOfX());
    trajectories.push_back(airResult.getYOfX());

    std::vector<std::string> labels{"Vacuum (DP5)", "Air Resistance (DP5)"};

    Visualizer::VisualizeMultiRealFunction(trajectories,
        "Projectile: Vacuum vs Air (v=" + std::to_string((int)velocity) + " m/s, 45 deg)",
        labels, 0, 1.05 * exactRange, 500, "projectile_vacuum_vs_air.mml");

    std::cout << ">> Scenario 1 complete - visualization saved\n";
}


/******************************************************************************
 * SCENARIO 2: Optimal Launch Angle
 * 
 * In vacuum, 45° gives maximum range. With air resistance?
 * The optimal angle is LOWER (typically 35-40°)!
 * Now computed with adaptive event detection - no tMax guessing needed.
 *****************************************************************************/
void Demo_OptimalAngle()
{
    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  SCENARIO 2: Optimal Launch Angle (Adaptive + Event Detection)   ║\n";
    std::cout << "╠═══════════════════════════════════════════════════════════════════╣\n";
    std::cout << "║  Testing 20-70 deg with vacuum and air resistance models         ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n\n";

    Real velocity = 50.0;
    Real initHeight = 0.0;
    std::cout << "Launch: v0 = " << velocity << " m/s, h0 = " << initHeight << " m\n\n";

    // Test angles 20° to 70° in 5° steps
    Vector<Real> angles(11);
    for (int i = 0; i < 11; i++)
        angles[i] = Utils::DegToRad(20.0 + i * 5.0);

    // ── Vacuum (analytical + numerical for verification) ──
    Projectile2DInVacuumODE vacuumSys;
    ProjectileMotionSolver2D vacuumSolver(vacuumSys);
    Vector<ProjectileTrajectory2D> vacuumResults = vacuumSolver.solveForAngles(angles, initHeight, velocity);

    // ── Air resistance ──
    Projectile2DWithAirResistanceODE airSys(0.01);
    ProjectileMotionSolver2D airSolver(airSys);
    Vector<ProjectileTrajectory2D> airResults = airSolver.solveForAngles(angles, initHeight, velocity);

    // Print comparison table
    std::cout << std::setw(10) << "Angle" 
              << std::setw(16) << "Vacuum (m)"
              << std::setw(16) << "Air (m)"
              << std::setw(14) << "Reduction" << "\n";
    std::cout << std::string(56, '-') << "\n";

    Real maxVacRange = 0, optVacAngle = 0;
    Real maxAirRange = 0, optAirAngle = 0;

    for (int i = 0; i < (int)angles.size(); i++)
    {
        Real deg = Utils::RadToDeg(angles[i]);
        Real vacRange = vacuumResults[i]._range;
        Real airRange = airResults[i]._range;
        Real reduction = (1 - airRange / vacRange) * 100;

        std::cout << std::fixed << std::setprecision(1)
                  << std::setw(9) << deg << " deg"
                  << std::setprecision(2)
                  << std::setw(16) << vacRange
                  << std::setw(16) << airRange
                  << std::setprecision(1)
                  << std::setw(12) << reduction << "%\n";

        if (vacRange > maxVacRange) { maxVacRange = vacRange; optVacAngle = deg; }
        if (airRange > maxAirRange) { maxAirRange = airRange; optAirAngle = deg; }
    }

    std::cout << "\n  Vacuum optimal: " << optVacAngle << " deg -> " 
              << std::setprecision(2) << maxVacRange << " m\n";
    std::cout << "  Air optimal:    " << optAirAngle << " deg -> " 
              << maxAirRange << " m\n";
    std::cout << "  Angle shift:    " << (optVacAngle - optAirAngle) << " deg lower with drag!\n\n";

    // Verify vacuum against analytical
    Real exactRange45 = vacuumSys.CalcRange(Utils::DegToRad(45.0), initHeight, velocity);
    Real numRange45 = vacuumResults[5]._range;  // index 5 = 45°
    std::cout << "  Vacuum 45 deg verification: exact=" << std::setprecision(6) << exactRange45
              << ", numerical=" << numRange45
              << ", err=" << std::scientific << std::setprecision(2) 
              << std::abs(numRange45 - exactRange45) << "\n\n";

    // ── Visualize air trajectories ──
    std::vector<LinearInterpRealFunc> trajectories;
    std::vector<std::string> labels;
    for (auto& result : airResults)
    {
        trajectories.push_back(result.getYOfX());
        labels.push_back(std::to_string((int)Utils::RadToDeg(result._angle)) + " deg");
    }

    Visualizer::VisualizeMultiRealFunction(trajectories,
        "Optimal Angle with Air Resistance (v=" + std::to_string((int)velocity) + " m/s)",
        labels, 0, 1.1 * maxAirRange, 500, "projectile_angles.mml");

    std::cout << ">> Scenario 2 complete - visualization saved\n";
}


/******************************************************************************
 * SCENARIO 3: Baseball Surface Comparison
 * 
 * Different ball surfaces = different drag coefficients!
 * Using speed-dependent drag models with adaptive event detection.
 *****************************************************************************/
void Demo_BaseballTypes()
{
    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  SCENARIO 3: Baseball Surface Comparison                         ║\n";
    std::cout << "╠═══════════════════════════════════════════════════════════════════╣\n";
    std::cout << "║  Speed-dependent drag + event-detected ground impact             ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n\n";

    Real velocity = 45.0;       // ~100 mph exit velocity
    Real angle = Utils::DegToRad(35.0);
    Real initHeight = 1.0;      // bat height

    std::cout << "Baseball hit: v0 = " << velocity << " m/s (~100 mph), angle = 35 deg, h0 = "
              << initHeight << " m\n\n";

    // Three ball types with speed-dependent drag
    BaseballWithDragCoeffDependentOnSpeedODE smoothBall(BaseballType::Smooth);
    BaseballWithDragCoeffDependentOnSpeedODE normalBall(BaseballType::Normal);
    BaseballWithDragCoeffDependentOnSpeedODE roughBall(BaseballType::Rough);

    ProjectileMotionSolver2D smoothSolver(smoothBall);
    ProjectileMotionSolver2D normalSolver(normalBall);
    ProjectileMotionSolver2D roughSolver(roughBall);

    // Vacuum reference
    Projectile2DInVacuumODE vacuumSys;
    ProjectileMotionSolver2D vacuumSolver(vacuumSys);
    ProjectileTrajectory2D vacResult = vacuumSolver.solveAdaptive(angle, initHeight, velocity);

    // Adaptive solve with event detection
    ProjectileTrajectory2D smoothResult = smoothSolver.solveAdaptive(angle, initHeight, velocity);
    ProjectileTrajectory2D normalResult = normalSolver.solveAdaptive(angle, initHeight, velocity);
    ProjectileTrajectory2D roughResult  = roughSolver.solveAdaptive(angle, initHeight, velocity);

    std::cout << std::setw(12) << "Type" 
              << std::setw(14) << "Range (m)"
              << std::setw(14) << "Flight (s)"
              << std::setw(14) << "Max H (m)" << "\n";
    std::cout << std::string(54, '-') << "\n";

    auto printRow = [](const std::string& name, const ProjectileTrajectory2D& r) {
        // Find max height from trajectory
        Real maxH = 0;
        for (int i = 0; i < (int)r._yValues.size(); i++)
            if (r._yValues[i] > maxH) maxH = r._yValues[i];

        std::cout << std::fixed << std::setprecision(2)
                  << std::setw(12) << name
                  << std::setw(14) << r._range
                  << std::setw(14) << r._timeOfFlight
                  << std::setw(14) << maxH << "\n";
    };

    printRow("Vacuum", vacResult);
    printRow("Smooth", smoothResult);
    printRow("Normal", normalResult);
    printRow("Rough",  roughResult);

    std::cout << "\n  Air resistance reduces range by "
              << std::setprecision(1) << (1 - normalResult._range / vacResult._range) * 100
              << "% for a normal baseball!\n\n";

    // ── Visualize ──
    std::vector<LinearInterpRealFunc> trajectories;
    trajectories.push_back(vacResult.getYOfX());
    trajectories.push_back(smoothResult.getYOfX());
    trajectories.push_back(normalResult.getYOfX());
    trajectories.push_back(roughResult.getYOfX());

    std::vector<std::string> labels{"Vacuum", "Smooth ball", "Normal ball", "Rough ball"};

    Visualizer::VisualizeMultiRealFunction(trajectories,
        "Baseball Trajectories by Surface Type",
        labels, 0, 1.1 * vacResult._range, 500, "projectile_baseball_types.mml");

    std::cout << ">> Scenario 3 complete - visualization saved\n";
}


/******************************************************************************
 * SCENARIO 4: Altitude-Dependent Air Density
 * 
 * For long-range artillery: air density decreases with altitude.
 * Comparing isothermal vs adiabatic (standard atmosphere) models.
 *****************************************************************************/
void Demo_AltitudeDrag()
{
    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  SCENARIO 4: Altitude-Dependent Air Density                      ║\n";
    std::cout << "╠═══════════════════════════════════════════════════════════════════╣\n";
    std::cout << "║  Long-range artillery: isothermal vs adiabatic atmosphere         ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n\n";

    Real velocity = 300.0;      // m/s - artillery speed
    Real angle = Utils::DegToRad(45.0);
    Real initHeight = 0.0;

    std::cout << "Artillery launch: v0 = " << velocity << " m/s, angle = 45 deg\n\n";

    // Models
    Projectile2DInVacuumODE vacuumSys;
    Projectile2DWithAirResistanceODE constAirSys(0.001);
    Projectile2DChangingAirDensityODE isoSys(0.001);
    Projectile2DChangingAirDensityODE adiaSys(0.001, AirDensityModel::Adiabatic);

    ProjectileMotionSolver2D vacSolver(vacuumSys);
    ProjectileMotionSolver2D constSolver(constAirSys);
    ProjectileMotionSolver2D isoSolver(isoSys);
    ProjectileMotionSolver2D adiaSolver(adiaSys);

    ProjectileTrajectory2D vacResult   = vacSolver.solveAdaptive(angle, initHeight, velocity, 0.05);
    ProjectileTrajectory2D constResult = constSolver.solveAdaptive(angle, initHeight, velocity, 0.05);
    ProjectileTrajectory2D isoResult   = isoSolver.solveAdaptive(angle, initHeight, velocity, 0.05);
    ProjectileTrajectory2D adiaResult  = adiaSolver.solveAdaptive(angle, initHeight, velocity, 0.05);

    std::cout << std::setw(20) << "Model"
              << std::setw(14) << "Range (m)"
              << std::setw(14) << "Flight (s)" << "\n";
    std::cout << std::string(48, '-') << "\n";

    auto printRow = [](const std::string& name, const ProjectileTrajectory2D& r) {
        std::cout << std::fixed << std::setprecision(2)
                  << std::setw(20) << name
                  << std::setw(14) << r._range
                  << std::setw(14) << r._timeOfFlight << "\n";
    };

    printRow("Vacuum", vacResult);
    printRow("Constant density", constResult);
    printRow("Isothermal atm.", isoResult);
    printRow("Adiabatic (std.)", adiaResult);

    std::cout << "\n  At high altitude, thinner air means LONGER range!\n";
    std::cout << "  Isothermal range increase over constant: +"
              << std::setprecision(1) 
              << (isoResult._range / constResult._range - 1) * 100 << "%\n\n";

    // ── Visualize ──
    std::vector<LinearInterpRealFunc> trajectories;
    trajectories.push_back(vacResult.getYOfX());
    trajectories.push_back(constResult.getYOfX());
    trajectories.push_back(isoResult.getYOfX());
    trajectories.push_back(adiaResult.getYOfX());

    std::vector<std::string> labels{"Vacuum", "Constant air", "Isothermal", "Adiabatic (std)"};

    Visualizer::VisualizeMultiRealFunction(trajectories,
        "Artillery: Air Density Models at 300 m/s",
        labels, 0, 1.05 * vacResult._range, 500, "projectile_altitude_drag.mml");

    std::cout << ">> Scenario 4 complete - visualization saved\n";
}


/******************************************************************************
 * MAIN
 *****************************************************************************/
int main()
{
    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║              MML PROJECTILE LAUNCH SIMULATOR                     ║\n";
    std::cout << "╠═══════════════════════════════════════════════════════════════════╣\n";
    std::cout << "║                                                                   ║\n";
    std::cout << "║  Solver: DormandPrince5 (adaptive, 5th order)                    ║\n";
    std::cout << "║  Ground: IODESystemWithEvents (zero-crossing, ~1e-12)            ║\n";
    std::cout << "║                                                                   ║\n";
    std::cout << "║  Scenarios:                                                      ║\n";
    std::cout << "║    1. Vacuum vs Air + Euler vs DP5 precision showdown            ║\n";
    std::cout << "║    2. Optimal launch angle with drag                             ║\n";
    std::cout << "║    3. Baseball surface types (speed-dependent drag)              ║\n";
    std::cout << "║    4. Altitude-dependent air density (artillery)                 ║\n";
    std::cout << "║                                                                   ║\n";
    std::cout << "║  Physics: F = mg - kv^2 (quadratic drag)                         ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n";

    Demo_VacuumVsAir();
    Demo_OptimalAngle();
    Demo_BaseballTypes();
    Demo_AltitudeDrag();

    std::cout << "\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  ALL PROJECTILE SIMULATIONS COMPLETE                             ║\n";
    std::cout << "╠═══════════════════════════════════════════════════════════════════╣\n";
    std::cout << "║  Output files (.mml):                                            ║\n";
    std::cout << "║    projectile_vacuum_vs_air.mml                                  ║\n";
    std::cout << "║    projectile_angles.mml                                         ║\n";
    std::cout << "║    projectile_baseball_types.mml                                 ║\n";
    std::cout << "║    projectile_altitude_drag.mml                                  ║\n";
    std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n\n";

    return 0;
}
