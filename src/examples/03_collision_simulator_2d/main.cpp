/******************************************************************************
 * MML Example: 2D Collision Simulator - Kinetic Theory & Shock Waves
 * ============================================================================
 * 
 * Demonstrates MML's particle collision physics with advanced features:
 * 
 *   1. TWO-COLOR MIXING - 500 balls (blue/red) mixing from separate halves
 *   2. SHOCK WAVE       - 30,000+ balls with energetic core explosion
 * 
 * Physics: Elastic collisions conserve both momentum and kinetic energy
 *   - Exact collision time calculation for sub-timestep accuracy
 *   - Space subdivision for O(N) average collision detection
 *   - Parallel execution for large-scale simulations
 * 
 * Features:
 *   - Self-contained CollisionSimulator2D with spatial partitioning
 *   - Multi-threaded collision detection for 30,000+ particles
 *   - ParticleVisualizer2D integration for real-time animation
 * 
 * Build: cmake --build build --target Example03_CollisionSim2D
 * Run:   ./build/src/examples/Debug/Example03_CollisionSim2D
 * 
 *****************************************************************************/

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "mml/tools/Serializer.h"
#include "mml/tools/Visualizer.h"
#endif

// Self-contained collision physics (no MPL dependency)
#include "CollisionSimulator2D.h"

#include <iostream>
#include <iomanip>
#include <chrono>

using namespace MML;
using namespace Collision2D;


/******************************************************************************
 * SCENARIO 1: Two-Color Ball Mixing
 * 
 * 500 balls (250 blue on left, 250 red on right) with random velocities.
 * Watch them mix over 100 time steps - demonstrates diffusion!
 *****************************************************************************/
void Demo_TwoColorMixing()
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  SCENARIO 1: Two-Color Ball Mixing (500 balls)\n";
    std::cout << std::string(70, '=') << "\n\n";

    // Configuration
    const double width = 1000.0;
    const double height = 800.0;
    const int nBalls = 250;          // 250 per color = 500 total
    const double mass = 1.0;
    const double radius = 5.0;
    const double velocity = 50.0;    // Initial speed
    const double dt = 0.1;
    const int numSteps = 100;

    std::cout << "Configuration:\n";
    std::cout << "  Box: " << width << " x " << height << "\n";
    std::cout << "  Balls: " << 2*nBalls << " total (250 blue + 250 red)\n";
    std::cout << "  Radius: " << radius << ", Mass: " << mass << "\n";
    std::cout << "  Initial velocity: " << velocity << " m/s\n";
    std::cout << "  Time step: " << dt << " s, Steps: " << numSteps << "\n\n";

    // Create container with two halves
    BoxContainer2D box = ContainerFactory2D::CreateTwoHalves(
        width, height,
        "Blue", nBalls, mass, radius, velocity,  // Left half
        "Red",  nBalls, mass, radius, velocity   // Right half
    );

    std::cout << "Created " << box.NumBalls() << " balls\n";

    // Create simulator with space subdivision (10x10 grid)
    // Grid should be larger than typical ball diameter
    int gridSize = 10;
    CollisionSimulator2D simulator(box, gridSize, gridSize);

    std::cout << "Running simulation with space subdivision (" 
              << gridSize << "x" << gridSize << " grid)...\n";

    auto start = std::chrono::high_resolution_clock::now();

    // Run simulation (use Fast mode for 500 balls - multithread overhead not worth it)
    SimResultsCollSim2D results = simulator.Simulate(numSteps, dt, RunTypeFast, true);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "\nSimulation completed in " << duration.count() << " ms\n";

    // Statistics
    double minSpeed, maxSpeed, avgSpeed, speedDev;
    results.CalcAllBallsStatistic(numSteps - 1, minSpeed, maxSpeed, avgSpeed, speedDev);
    
    std::cout << "\nFinal Statistics:\n";
    std::cout << "  Min speed: " << std::fixed << std::setprecision(2) << minSpeed << "\n";
    std::cout << "  Max speed: " << maxSpeed << "\n";
    std::cout << "  Avg speed: " << avgSpeed << " (σ = " << speedDev << ")\n";

    // Save and visualize
    std::string filename = "results/collision_two_color_mixing.txt";
    bool success = simulator.Serialize(filename, results, dt, 1);

    if (success)
    {
        std::cout << "\nSimulation saved to: " << filename << "\n";
        std::cout << "Launching ParticleVisualizer2D...\n";
        Visualizer::VisualizeParticleSimulation2D("collision_two_color_mixing.txt");
    }
    else
    {
        std::cout << "\nFailed to save simulation!\n";
    }

    std::cout << "\n✓ Two-color mixing complete!\n";
}


/******************************************************************************
 * SCENARIO 2: Shock Wave Simulation
 * 
 * 30,000 balls with an energetic core at center.
 * Requires parallel execution and space subdivision for performance.
 * 
 * Physics: Hot core expands into cold gas, creating shock wave.
 * The wave propagates outward as energy transfers through collisions.
 *****************************************************************************/
void Demo_ShockWave()
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  SCENARIO 2: Shock Wave Simulation (30,000 balls)\n";
    std::cout << std::string(70, '=') << "\n\n";

    // Configuration
    const double width = 2000.0;
    const double height = 1600.0;
    
    const int numSurrounding = 29900;   // Cold gas
    const double mass1 = 1.0;
    const double rad1 = 2.0;
    const double vel1 = 5.0;            // Low speed (cold)

    const int numEnergetic = 100;       // Hot core
    const double coreRadius = 50.0;
    const double mass2 = 1.0;
    const double rad2 = 2.0;
    const double vel2 = 500.0;          // High speed (hot!)

    const double dt = 0.01;
    const int numSteps = 200;

    std::cout << "Configuration:\n";
    std::cout << "  Box: " << width << " x " << height << "\n";
    std::cout << "  Surrounding gas: " << numSurrounding << " balls (v=" << vel1 << ")\n";
    std::cout << "  Energetic core: " << numEnergetic << " balls (v=" << vel2 << ")\n";
    std::cout << "  Core radius: " << coreRadius << "\n";
    std::cout << "  Time step: " << dt << " s, Steps: " << numSteps << "\n\n";

    std::cout << "Creating shock wave configuration...\n";

    BoxContainer2D box = ContainerFactory2D::CreateShockWave(
        width, height,
        numSurrounding, mass1, rad1, vel1,
        numEnergetic, coreRadius, mass2, rad2, vel2
    );

    std::cout << "Created " << box.NumBalls() << " balls\n";

    // Create simulator with fine grid for 30K balls
    // Grid cells should be ~4x ball diameter for good performance
    int gridRows = static_cast<int>(height / (4 * rad1));
    int gridCols = static_cast<int>(width / (4 * rad1));
    
    // Cap at reasonable size
    gridRows = std::min(gridRows, 200);
    gridCols = std::min(gridCols, 250);

    std::cout << "Using " << gridRows << "x" << gridCols << " spatial grid\n";

    CollisionSimulator2D simulator(box, gridRows, gridCols);

    std::cout << "Running PARALLEL simulation...\n";
    std::cout << "  (This will use all CPU cores for collision detection)\n\n";

    auto start = std::chrono::high_resolution_clock::now();

    // Run with multi-threading for 30K balls
    SimResultsCollSim2D results = simulator.Simulate(
        numSteps, dt, RunTypeFastMultithread, true);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    std::cout << "\nSimulation completed in " << duration.count() << " seconds\n";

    // Statistics at different times
    std::cout << "\nEnergy propagation (avg speed evolution):\n";
    std::cout << "  Step    Avg Speed\n";
    std::cout << "  ----    ---------\n";
    
    for (int step = 0; step < numSteps; step += 40)
    {
        double minSpeed, maxSpeed, avgSpeed, speedDev;
        results.CalcAllBallsStatistic(step, minSpeed, maxSpeed, avgSpeed, speedDev);
        std::cout << "  " << std::setw(4) << step << "    " 
                  << std::fixed << std::setprecision(2) << avgSpeed << "\n";
    }

    // Save and visualize
    std::string filename = "results/collision_shock_wave.txt";
    
    // Save every 2nd frame to reduce file size
    bool success = simulator.Serialize(filename, results, dt, 2);

    if (success)
    {
        std::cout << "\nSimulation saved to: " << filename << "\n";
        std::cout << "Launching ParticleVisualizer2D...\n";
        Visualizer::VisualizeParticleSimulation2D("collision_shock_wave.txt");
    }
    else
    {
        std::cout << "\nFailed to save simulation!\n";
    }

    std::cout << "\n✓ Shock wave simulation complete!\n";
}


/******************************************************************************
 * MAIN
 *****************************************************************************/
int main()
{
    std::cout << "\n";
    std::cout << "======================================================================\n";
    std::cout << "         MML 2D COLLISION SIMULATOR - ADVANCED PHYSICS\n";
    std::cout << "======================================================================\n\n";
    
    std::cout << "   Self-contained collision physics engine demonstrating:\n";
    std::cout << "   - Exact elastic collision calculation\n";
    std::cout << "   - Spatial subdivision for O(N) performance\n";
    std::cout << "   - Parallel execution for large simulations\n";
    std::cout << "   - Real-time visualization with ParticleVisualizer2D\n";
    std::cout << "======================================================================\n";

    // Run scenarios
    Demo_TwoColorMixing();
    Demo_ShockWave();

    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  ALL COLLISION SIMULATIONS COMPLETE!\n";
    std::cout << std::string(70, '=') << "\n\n";

    std::cout << "Output files in results/:\n";
    std::cout << "  - collision_two_color_mixing.txt (500 balls, 100 steps)\n";
    std::cout << "  - collision_shock_wave.txt (30,000 balls, 200 steps)\n\n";

    std::cout << "Visualization:\n";
    std::cout << "  ParticleVisualizer2D launched automatically for each simulation.\n";
    std::cout << "  Use Play/Pause and speed controls to explore the physics!\n\n";

    return 0;
}
