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

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/tools/Serializer.h"
#include "mml/tools/Visualizer.h"
#include "mml/tools/ConsolePrinter.h"

#include "mpl/Base/SolarSystem.h"
#include "mpl/Gravity/NBodySimulator.h"
#endif

#include <iostream>
#include <iomanip>
#include <fstream>

using namespace MML;
using namespace MPL;

void Demo_StarClusterCollision()
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  SCENARIO: COLLISION OF STAR CLUSTERS\n";
    std::cout << std::string(70, '=') << "\n\n";

    std::cout << "Two 100-body star clusters swing past each other!\n";
    std::cout << "  Cluster A (Cyan):    100 stars + central mass, moving RIGHT\n";
    std::cout << "  Cluster B (Magenta): 100 stars + central mass, moving LEFT\n";
    std::cout << "  Initial separation: 800 units\n";
    std::cout << "  Impact parameter: 250 units (hyperbolic flyby)\n";
    std::cout << "  Approach speed: 5.0 units/time (relative)\n\n";

    // Create star cluster collision configuration
    NBodyGravitySimConfig config = NBodyGravityConfigGenerator::Config4_StarClusterCollision(
        100,      // bodies per cluster
        800.0,    // separation between cluster centers
        250.0,    // impact parameter - offset for flyby
        5.0,      // approach speed (doubled for faster encounter!)
        1.0       // G
    );

    std::cout << "Configuration Created:\n";
    std::cout << "  Total bodies: " << config.NumBodies() << " (2 clusters x 101)\n";
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
    std::cout << "  Complexity: O(N^2) = " << config.NumBodies() * config.NumBodies() 
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
    std::string outputFile = "results/star_cluster_collision.mml";
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

    std::cout << "\n  Star Cluster Collision simulation complete!\n";
}