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

void Demo_Solar_system()
{
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "  SCENARIO: Solar System Simulation\n";
    std::cout << std::string(70, '=') << "\n\n";

    // Create solar system configuration (uses real astronomical data!)
    NBodyGravitySimConfig config = NBodyGravityConfigGenerator::Config2_Solar_system();

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
    std::string outputFile = "results/solar_system_trajectories.mml";
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

    std::cout << "\n  Solar system simulation complete!\n";
}

