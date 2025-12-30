/**
 * @file show_particle_visualizer_3d.cpp
 * @brief Demonstrates 3D particle simulation visualization
 * 
 * Uses the MML Visualizer infrastructure to display 3D particle simulations.
 * Loads pre-generated simulation data from the visualization_data subfolder.
 * Cross-platform: works on Windows, Linux, and macOS.
 */

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "tools/Visualizer.h"
#endif

#include <filesystem>

using namespace MML;

void Show_Particle_Visualizer_3D_Examples()
{
    std::cout << "\n=== 3D Particle Simulation Visualization Examples ===\n\n";

    // Get the path to the visualization data folder
    std::string dataFolder = "src/visualization_examples/visualization_data/";
    std::string dataFile = "particle_visualizer_3d_1.txt";
    
    // First, try to find the file relative to current working directory
    std::filesystem::path filePath = dataFolder + dataFile;
    
    // If not found, try looking in common alternative locations
    if (!std::filesystem::exists(filePath)) {
        std::vector<std::string> searchPaths = {
            "../" + dataFolder + dataFile,
            "../../" + dataFolder + dataFile,
            "../../../" + dataFolder + dataFile,
            dataFile
        };
        
        for (const auto& path : searchPaths) {
            if (std::filesystem::exists(path)) {
                filePath = path;
                break;
            }
        }
    }
    
    if (std::filesystem::exists(filePath)) {
        std::cout << "Loading 3D particle simulation from: " << filePath << "\n";
        std::cout << "This shows a solar system simulation with 9 bodies.\n\n";
        
        // Copy the file to results folder where the visualizer expects it
        std::filesystem::path resultsPath = GetResultFilesPath() + dataFile;
        
        try {
            std::filesystem::copy_file(filePath, resultsPath, 
                                       std::filesystem::copy_options::overwrite_existing);
            
            // Now visualize using the Visualizer class
            Visualizer::VisualizeParticleSimulation3D(dataFile);
            
        } catch (const std::exception& e) {
            std::cerr << "Error copying file: " << e.what() << "\n";
            std::cerr << "Attempting direct visualization...\n";
            
            std::cout << "Note: Direct visualization from source location not implemented.\n";
            std::cout << "Please manually copy the file to: " << GetResultFilesPath() << "\n";
        }
    } else {
        std::cout << "3D particle data file not found.\n";
        std::cout << "Expected location: " << filePath << "\n";
        std::cout << "Searched in: current directory and parent directories\n\n";
        
        std::cout << "To generate sample data, run a gravity simulation:\n";
        std::cout << "  - See NBodySimulator.h for N-body gravity simulations\n";
        std::cout << "  - Or collision_simulator_3d.cpp for 3D collision examples\n";
        
        // Show expected file format
        std::cout << "\nExpected file format:\n";
        std::cout << "  PARTICLE_SIMULATION_DATA_3D\n";
        std::cout << "  Width: <W>\n";
        std::cout << "  Height: <H>\n";
        std::cout << "  Depth: <D>\n";
        std::cout << "  NumBalls: <N>\n";
        std::cout << "  Ball_1 <Color> <Radius>\n";
        std::cout << "  Ball_2 <Color> <Radius>\n";
        std::cout << "  ...\n";
        std::cout << "  NumSteps: <M>\n";
        std::cout << "  Step <i> <time>\n";
        std::cout << "  <id> <x> <y> <z>\n";
        std::cout << "  ...\n";
    }

    std::cout << "\n3D particle visualization example complete!\n";
}
