///////////////////////////////////////////////////////////////////////////////////////////
// Chapter 02 - Particle Visualization Demos
// 
// This file demonstrates MML's particle simulation visualization capabilities.
// Particle visualizers animate multiple particles over time (2D and 3D).
//
// Demos use pre-generated simulation data files.
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/tools/Visualizer.h"
#endif

#include <filesystem>

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
//                           PARTICLE VISUALIZATION DEMOS
///////////////////////////////////////////////////////////////////////////////////////////

void Chapter02_Demo_Particles()
{
	std::cout << "\n=== Chapter 02: Particle Visualization ===" << std::endl;

	// Get the path to data files (relative to source file location)
	std::filesystem::path dataPath = std::filesystem::path(__FILE__).parent_path() / "data";

	// Demo 1: 2D Particle Simulation
	// 100 red balls bouncing in a 1000x1000 box
	{
		std::cout << "  Visualizing 2D Particle Simulation (100 balls)..." << std::endl;
		std::filesystem::path file2D = dataPath / "particle_visualizer_2d_1.mml";
		
		if (std::filesystem::exists(file2D)) {
			auto result = Visualizer::VisualizeParticleSimulation2D(file2D.string());
			if (!result.success) {
				std::cerr << "    Failed: " << result.errorMessage << std::endl;
			}
		} else {
			std::cerr << "    Data file not found: " << file2D << std::endl;
		}
	}

	// Demo 2: 3D Particle Simulation
	{
		std::cout << "  Visualizing 3D Particle Simulation..." << std::endl;
		std::filesystem::path file3D = dataPath / "particle_visualizer_3d_1.mml";
		
		if (std::filesystem::exists(file3D)) {
			auto result = Visualizer::VisualizeParticleSimulation3D(file3D.string());
			if (!result.success) {
				std::cerr << "    Failed: " << result.errorMessage << std::endl;
			}
		} else {
			std::cerr << "    Data file not found: " << file3D << std::endl;
		}
	}

	std::cout << "  Particle demos complete!" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
//                           MAIN FUNCTION FOR THIS MODULE
///////////////////////////////////////////////////////////////////////////////////////////

void Chapter02_ParticleVisualization()
{
	std::cout << "\n========================================" << std::endl;
	std::cout << "Chapter 02: Particle Visualization" << std::endl;
	std::cout << "========================================" << std::endl;

	Chapter02_Demo_Particles();

	std::cout << "\nAll particle visualizations complete!" << std::endl;
}
