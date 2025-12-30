/**
 * @file show_particle_visualizer_2d.cpp
 * @brief Demonstrates 2D particle simulation visualization
 * 
 * Uses the MML Visualizer infrastructure to display 2D particle simulations.
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

void Show_Particle_Visualizer_2D_Examples()
{
    std::cout << "\n=== 2D Particle Simulation Visualization Examples ===\n\n";

    // Get the path to the visualization data folder
    // The data files are located relative to the source directory
    // We need to construct the path based on where the executable is run from
    
    std::string dataFolder = "src/visualization_examples/visualization_data/";
    std::string dataFile = "particle_visualizer_2d_1.txt";
    
    // First, try to find the file relative to current working directory
    std::filesystem::path filePath = dataFolder + dataFile;
    
    // If not found, try looking in common alternative locations
    if (!std::filesystem::exists(filePath)) {
        // Try relative to parent directories (common when running from build folder)
        std::vector<std::string> searchPaths = {
            "../" + dataFolder + dataFile,
            "../../" + dataFolder + dataFile,
            "../../../" + dataFolder + dataFile,
            // Also try if we're in the results folder already
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
        std::cout << "Loading 2D particle simulation from: " << filePath << "\n";
        std::cout << "This shows a 2D collision simulation with 100 particles.\n\n";
        
        // Copy the file to results folder where the visualizer expects it
        std::filesystem::path resultsPath = GetResultFilesPath() + dataFile;
        
        try {
            std::filesystem::copy_file(filePath, resultsPath, 
                                       std::filesystem::copy_options::overwrite_existing);
            
            // Now visualize using the Visualizer class
            Visualizer::VisualizeParticleSimulation2D(dataFile);
            
        } catch (const std::exception& e) {
            std::cerr << "Error copying file: " << e.what() << "\n";
            std::cerr << "Attempting direct visualization...\n";
            
            // Try to visualize directly from the source location
            // This requires the visualizer to accept full paths
            std::cout << "Note: Direct visualization from source location not implemented.\n";
            std::cout << "Please manually copy the file to: " << GetResultFilesPath() << "\n";
        }
    } else {
        std::cout << "2D particle data file not found.\n";
        std::cout << "Expected location: " << filePath << "\n";
        std::cout << "Searched in: current directory and parent directories\n\n";
        
        std::cout << "To generate sample data, run a collision simulation:\n";
        std::cout << "  - See Chapter 4 (collision_simulator_2d.cpp) for examples\n";
        std::cout << "  - Or create your own PARTICLE_SIMULATION_DATA_2D file\n";
        
        // Show expected file format
        std::cout << "\nExpected file format:\n";
        std::cout << "  PARTICLE_SIMULATION_DATA_2D\n";
        std::cout << "  NumBalls: <N>\n";
        std::cout << "  Ball_1 <Color> <Radius>\n";
        std::cout << "  Ball_2 <Color> <Radius>\n";
        std::cout << "  ...\n";
        std::cout << "  NumSteps: <M>\n";
        std::cout << "  Step <i> <time>\n";
        std::cout << "  <x1> <y1>\n";
        std::cout << "  <x2> <y2>\n";
        std::cout << "  ...\n";
    }

    std::cout << "\n2D particle visualization example complete!\n";
}
