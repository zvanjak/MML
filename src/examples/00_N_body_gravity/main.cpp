/******************************************************************************
 * MML PRIME EXAMPLE: N-Body Gravity Simulator
 * ============================================================================
 * 
 * This is the flagship demonstration of MML's numerical computing power!
 * 
 * Demonstrates MML's gravitational N-body simulation capabilities with three
 * impressive scenarios:
 * 
 *   1. SOLAR SYSTEM          - Our actual solar system with real masses
 *   2. MANY BODIES           - 100 bodies orbiting a massive central object  
 *   3. STAR CLUSTER COLLISION - TWO 100-body clusters colliding! (SPECTACULAR!)
 * 
 * Physics: Newton's law of universal gravitation
 *   F = G * m1 * m2 / r^2  (attractive force along line joining masses)
 * 
 * Numerical Methods:
 *   - Euler method (simple, illustrative)
 *   - RK5 Cash-Karp (adaptive, production-quality)
 * 
 * Output: Trajectory files for visualization with MML Visualizer
 * 
 * Build: cmake --build build --target Example00_NBodyGravity
 * Run:   ./build/src/examples/Debug/Example00_NBodyGravity
 * 
 *****************************************************************************/

#include "MMLBase.h"
#include "mml/tools/Serializer.h"
#include "mml/tools/Visualizer.h"
#include "mml/tools/ConsolePrinter.h"

// Self-contained N-body simulation (no MPL dependency)
#include "NBodyGravity.h"

#include <iostream>
#include <iomanip>
#include <fstream>

using namespace MML;
using namespace NBody;

/******************************************************************************
 * SCENARIO 1: Solar System Simulation
 * 
 * Simulates the actual solar system with:
 *   - Sun (center)
 *   - All 8 planets with real masses and orbital distances
 *   - Circular initial orbits
 * 
 * Units: million km, Jupiter masses, years
 *****************************************************************************/
void Demo_SolarSystem()
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  SCENARIO 1: Solar System Simulation\n";
    std::cout << std::string(70, '=') << "\n\n";

    // Create solar system configuration (uses real astronomical data!)
    NBodyGravitySimConfig config = NBodyGravityConfigGenerator::Config1_Solar_system();

    std::cout << "Solar System Configuration:\n";
    std::cout << "  Bodies: " << config.NumBodies() << " (Sun + 8 planets)\n";
    std::cout << "  Gravitational constant G = " << std::scientific << config.G() << "\n";
    std::cout << "  Units: million km, Jupiter masses, years\n\n";

    // Create simulator
    NBodyGravitySimulator solver(config);

    // Simulation parameters
    Real duration = 2.0;        // 2 years (to see inner planet orbits clearly)
    const int steps = 1001;
    const Real dt = duration / steps;

    std::cout << "Simulation Parameters:\n";
    std::cout << "  Duration: " << duration << " years\n";
    std::cout << "  Time step: " << std::fixed << std::setprecision(4) << dt << " years\n";
    std::cout << "  Steps: " << steps << "\n\n";

    // Solve using Euler method
    std::cout << "Running Euler integration... " << std::flush;
    NBodyGravitySimulationResults result = solver.SolveEuler(dt, steps);
    std::cout << "Done!\n\n";

    // Output trajectory data
    std::string outputFile = "results/solar_system_trajectories.txt";
    std::cout << "Writing trajectories to: " << outputFile << "\n";

    std::ofstream ofs(outputFile);
    ofs << "# Solar System N-Body Simulation\n";
    ofs << "# Time(years) Sun_x Sun_y Mercury_x Mercury_y Venus_x Venus_y Earth_x Earth_y Mars_x Mars_y Jupiter_x Jupiter_y\n";
    ofs << "# Units: million km\n";
    ofs << std::fixed << std::setprecision(4);

    for (int i = 0; i < result.NumSteps(); i += 10)  // Every 10th step
    {
        ofs << result.Time(i);
        for (int body = 0; body < std::min(6, result.NumBodies()); body++)
        {
            ofs << " " << result.State(i).Pos(body).X() 
                << " " << result.State(i).Pos(body).Y();
        }
        ofs << "\n";
    }
    ofs.close();

    // Energy conservation check
    std::cout << "\nEnergy Conservation Check:\n";
    std::cout << "  Initial total energy: " << std::scientific << result.getInitialEnergy() << "\n";
    std::cout << "  Final total energy:   " << result.getFinalEnergy() << "\n";
    std::cout << "  Relative error: " << result.getRelativeEnergyError() * 100 << "%\n";

    // Visualization - launches WPF window, waits for user to close it
    result.VisualizeAsParamCurve("solar_system", Vector<int>{ 1, 2, 3, 4, 5, 6, 7, 8 });

    // Particle animation - 3D animated visualization
    result.VisualizeAsParticleSimulation("solar_system_anim", Vector<int>{ 0, 1, 2, 3, 4, 5, 6, 7, 8 }, dt);

    std::cout << "\n✓ Solar system simulation complete!\n";
}


/******************************************************************************
 * SCENARIO 2: Many Bodies Around Central Mass
 * 
 * 100 small bodies orbiting a massive central object.
 * Demonstrates handling of larger N-body systems.
 *****************************************************************************/
void Demo_ManyBodies()
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  SCENARIO 2: 100 Bodies Around Central Mass\n";
    std::cout << std::string(70, '=') << "\n\n";

    // Create configuration with many bodies
    NBodyGravitySimConfig config = NBodyGravityConfigGenerator::Config2_N_bodies_around_massive_object();

    std::cout << "Many-Body Configuration:\n";
    std::cout << "  Bodies: " << config.NumBodies() << "\n";
    std::cout << "  Central mass: " << config.Mass(0) << "\n";
    std::cout << "  Orbiting bodies: " << config.NumBodies() - 1 << "\n\n";

    // Create simulator
    NBodyGravitySimulator solver(config);

    // Shorter duration, fewer steps (computationally intensive!)
    Real duration = 50.0;
    const int steps = 501;
    const Real dt = duration / steps;

    std::cout << "Simulation Parameters:\n";
    std::cout << "  Duration: " << duration << " time units\n";
    std::cout << "  Time step: " << dt << "\n";
    std::cout << "  Steps: " << steps << "\n";
    std::cout << "  (Note: O(N²) complexity, " << config.NumBodies() * config.NumBodies() 
              << " force calculations per step)\n\n";

    // Solve
    std::cout << "Running Euler integration... " << std::flush;
    NBodyGravitySimulationResults result = solver.SolveEuler(dt, steps);
    std::cout << "Done!\n\n";

    // Output
    std::string outputFile = "results/many_body_trajectories.txt";
    std::cout << "Writing trajectories to: " << outputFile << "\n";

    std::ofstream ofs(outputFile);
    ofs << "# Many-Body Simulation (100+ bodies)\n";
    ofs << "# Time CentralMass_x CentralMass_y CentralMass_z (+ first 10 orbiting bodies)\n";
    ofs << std::fixed << std::setprecision(4);

    for (int i = 0; i < result.NumSteps(); i += 5)
    {
        ofs << result.Time(i);
        for (int body = 0; body < std::min(11, result.NumBodies()); body++)
        {
            ofs << " " << result.State(i).Pos(body).X() 
                << " " << result.State(i).Pos(body).Y()
                << " " << result.State(i).Pos(body).Z();
        }
        ofs << "\n";
    }
    ofs.close();

    // Visualization - show ALL bodies!
    Vector<int> bodiesToShow(config.NumBodies());
    for (int i = 0; i < config.NumBodies(); i++) bodiesToShow[i] = i;
    result.VisualizeAsParamCurve("many_bodies", bodiesToShow);

    // Particle animation - 3D animated visualization of all bodies
    result.VisualizeAsParticleSimulation("many_bodies_anim", bodiesToShow, dt);

    std::cout << "\n✓ Many-body simulation complete!\n";
}


/******************************************************************************
 * SCENARIO 3: COLLISION OF STAR CLUSTERS (THE SPECTACULAR ONE!)
 * 
 * Two separate 100-body star clusters approach and COLLIDE!
 * Each cluster has:
 *   - A central massive body (black hole / dense core)
 *   - 100 stars in bound orbits around it
 * 
 * The clusters start 800 units apart and move toward each other.
 * Watch gravitational dynamics unfold as they merge!
 * 
 * This is the most visually impressive scenario - demonstrates:
 *   - Large-scale gravitational interactions
 *   - Chaotic dynamics during collision
 *   - Conservation laws in complex systems
 *   - Beautiful trajectory visualization
 *****************************************************************************/
void Demo_StarClusterCollision()
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  SCENARIO 3: ★ COLLISION OF STAR CLUSTERS ★\n";
    std::cout << std::string(70, '=') << "\n\n";

    std::cout << "Two 100-body star clusters swing past each other!\n";
    std::cout << "  • Cluster A (Cyan):    100 stars + central mass, moving RIGHT\n";
    std::cout << "  • Cluster B (Magenta): 100 stars + central mass, moving LEFT\n";
    std::cout << "  • Initial separation: 800 units\n";
    std::cout << "  • Impact parameter: 250 units (hyperbolic flyby)\n";
    std::cout << "  • Approach speed: 5.0 units/time (relative)\n\n";

    // Create star cluster collision configuration
    NBodyGravitySimConfig config = NBodyGravityConfigGenerator::Config3_StarClusterCollision(
        100,      // bodies per cluster
        800.0,    // separation between cluster centers
        250.0,    // impact parameter - offset for flyby
        5.0,      // approach speed (doubled for faster encounter!)
        1.0       // G
    );

    std::cout << "Configuration Created:\n";
    std::cout << "  Total bodies: " << config.NumBodies() << " (2 clusters × 101)\n";
    std::cout << "  Central masses: 5000 each (white)\n";
    std::cout << "  Star masses: 1-4 (randomized)\n\n";

    // Create simulator
    NBodyGravitySimulator solver(config);

    // Simulation parameters - need long duration for approach + collision + aftermath
    Real duration = 600.0;     // Long enough for full collision dynamics
    const int steps = 3001;    // More steps for smooth animation
    const Real dt = duration / steps;

    std::cout << "Simulation Parameters:\n";
    std::cout << "  Duration: " << duration << " time units\n";
    std::cout << "  Time step: " << std::fixed << std::setprecision(4) << dt << "\n";
    std::cout << "  Steps: " << steps << "\n";
    std::cout << "  Complexity: O(N²) = " << config.NumBodies() * config.NumBodies() 
              << " force calculations per step\n\n";

    std::cout << "Running simulation... (this will take a moment)\n" << std::flush;
    
    // Use Euler for speed (RK5 would be more accurate but slower for 202 bodies)
    NBodyGravitySimulationResults result = solver.SolveEuler(dt, steps);
    std::cout << "Done! Simulated " << result.NumSteps() << " timesteps.\n\n";

    // Conservation laws verification
    std::cout << "Conservation Laws Check:\n";
    
    // Momentum (should be ~zero since clusters approach symmetrically)
    Vec3Cart p_initial = result.State(0).LinearMomentum();
    Vec3Cart p_final = result.State(result.NumSteps()-1).LinearMomentum();
    std::cout << "  Linear Momentum:\n";
    std::cout << "    Initial: (" << std::scientific << std::setprecision(3) 
              << p_initial.X() << ", " << p_initial.Y() << ", " << p_initial.Z() << ")\n";
    std::cout << "    Final:   (" << p_final.X() << ", " << p_final.Y() << ", " << p_final.Z() << ")\n";

    // Energy
    std::cout << "  Total Energy:\n";
    std::cout << "    Initial: " << std::scientific << result.getInitialEnergy() << "\n";
    std::cout << "    Final:   " << result.getFinalEnergy() << "\n";
    std::cout << "    Drift:   " << std::fixed << std::setprecision(2) 
              << result.getRelativeEnergyError() * 100 << "%\n";
    std::cout << "    (Note: Euler method has energy drift; RK5 would preserve better)\n\n";

    // Output trajectory file
    std::string outputFile = "results/star_cluster_collision.txt";
    std::cout << "Writing trajectories to: " << outputFile << "\n";

    std::ofstream ofs(outputFile);
    ofs << "# Star Cluster Collision - N-Body Simulation\n";
    ofs << "# " << config.NumBodies() << " bodies (2 clusters of ~101 each)\n";
    ofs << "# Time followed by x,y,z for first 20 bodies\n";
    ofs << std::fixed << std::setprecision(4);

    for (int i = 0; i < result.NumSteps(); i += 10)  // Every 10th step
    {
        ofs << result.Time(i);
        for (int body = 0; body < std::min(20, result.NumBodies()); body++)
        {
            ofs << " " << result.State(i).Pos(body).X() 
                << " " << result.State(i).Pos(body).Y()
                << " " << result.State(i).Pos(body).Z();
        }
        ofs << "\n";
    }
    ofs.close();

    // Visualization - show ALL bodies for the spectacular effect!
    std::cout << "\nLaunching visualization...\n";
    std::cout << "  (Prepare for a SPECTACULAR view of 202 trajectories!)\n\n";

    Vector<int> allBodies(config.NumBodies());
    for (int i = 0; i < config.NumBodies(); i++) allBodies[i] = i;

    // 3D trajectory curves (the "spaghetti" plot)
    result.VisualizeAsParamCurve("star_cluster_collision", allBodies);

    // Particle animation
    result.VisualizeAsParticleSimulation("star_cluster_anim", allBodies, dt);

    std::cout << "\n★ Star Cluster Collision simulation complete! ★\n";
}


/******************************************************************************
 * MAIN
 *****************************************************************************/
int main()
{
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║         MML N-BODY GRAVITATIONAL SIMULATOR                           ║\n";
    std::cout << "║         =====================================                        ║\n";
    std::cout << "║                                                                      ║\n";
    std::cout << "║   Simulating gravitational interactions between multiple bodies      ║\n";
    std::cout << "║   using Newton's law of universal gravitation.                       ║\n";
    std::cout << "║                                                                      ║\n";
    std::cout << "║   F = G * m₁ * m₂ / r²                                               ║\n";
    std::cout << "║                                                                      ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════════╝\n";

    // Run all scenarios
    Demo_SolarSystem();
    Demo_ManyBodies();
    Demo_StarClusterCollision();  // THE SPECTACULAR ONE!

    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  ★ ALL SIMULATIONS COMPLETE! ★\n";
    std::cout << std::string(70, '=') << "\n\n";

    std::cout << "Output files created in results/:\n";
    std::cout << "  - solar_system_trajectories.txt\n";
    std::cout << "  - many_body_trajectories.txt\n";
    std::cout << "  - star_cluster_collision.txt  ← THE SPECTACULAR ONE!\n\n";

    std::cout << "MML Visualization:\n";
    std::cout << "  Trajectories automatically visualized using MML's built-in\n";
    std::cout << "  Visualizer class with gnuplot integration.\n";
    std::cout << "  See: VisualizeAsParamCurve() and VisualizeAsParticleSimulation()\n\n";

    return 0;
}