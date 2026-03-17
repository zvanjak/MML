/******************************************************************************
 * SCENARIO 2: Many Bodies Around Central Mass
 * 
 * 100 small bodies orbiting a massive central object.
 * Demonstrates handling of larger N-body systems.
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

void Demo_ManyBodies()
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  SCENARIO: 100 Bodies Around Central Mass\n";
    std::cout << std::string(70, '=') << "\n\n";

    // Create configuration with many bodies
    NBodyGravitySimConfig config = NBodyGravityConfigGenerator::Config3_N_bodies_around_massive_object();

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
    std::cout << "  (Note: O(N^2) complexity, " << config.NumBodies() * config.NumBodies() 
              << " force calculations per step)\n\n";

    // Solve
    std::cout << "Running Euler integration... " << std::flush;
    NBodyGravitySimulationResults result = solver.SolveEuler(dt, steps);
    std::cout << "Done!\n\n";

    // Output
    std::string outputFile = "results/many_body_trajectories.mml";
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

    std::cout << "\n  Many-body simulation complete!\n";
}