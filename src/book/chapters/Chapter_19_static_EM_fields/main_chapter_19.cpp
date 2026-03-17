#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector/VectorTypes.h"
#include "base/Geometry/Geometry3DBodies.h"
#include "tools/Serializer.h"
#include "tools/Visualizer.h"

#include "mpl/Electromagnetism/GenericChargeSimulator.h"
#endif

#include <cmath>
#include <memory>

using namespace MML;
using namespace MPL;

// Verify Gauss' law for a point charge
// integracija po liniji, calculate work done by electric field

// Biot-Savart law - calculate exact magnetic field for a wire loop

// egzaktno magnetsko polje zavojnice, s N navoja, promjerom D, i strujom I
// Biot-Savartov zakon ... za svaku tocku integrirati po cijeloj liniji zavojnice

// TREBA MI VIZUALIZACIJA SILNICA b POLJA!!!
// dobro je sto su zatvorene petlje

// CHARGE DISTRIBUTION DEMOS
// First demo - fill given geometry randomly with electric charges
// and simulate their final distribution

// Second demo - fill given geometry with electric charges one by one
// and simulate their final distribution after each charge addition

// What happens if, after distributing all charges, we attach long wire to the geometry
// and connect it to the ground? Will the charges redistribute?

//////////////////////////////////////////////////////////////////////////////
// Demo functions - all use the same GenericChargeSimulator with different geometries!
//////////////////////////////////////////////////////////////////////////////

// Helper function to run simulation and visualize
void runChargeSimulationDemo(
	GenericChargeSimulator& sim,
	const std::string& geometryName,
	const std::string& fileName,
	const std::string& color,
	int numParticles,
	Real particleSizeScale = 1.0,
	int maxSteps = 2000)
{
	std::cout << "\n=== Demo: Charge Distribution in " << geometryName << " ===" << std::endl;
	std::cout << "Number of particles: " << numParticles << std::endl;

	// Simulation parameters
	sim.setTimeStep(0.1);
	sim.setDamping(0.98);
	sim.setCoulombConstant(10000.0);

	// Initialize with identical charges
	Real charge = 1.0;
	Real mass = 1.0;
	sim.initializeRandomParticles(numParticles, charge, mass);

	// Run simulation
	int saveEveryN = 1;
	auto history = sim.runToEquilibrium(maxSteps, 1e-12, saveEveryN, true);

	// Debug output
	std::cout << "\nDebug: Particle 0 positions:" << std::endl;
	for (size_t frame = 0; frame < std::min(history[0].size(), size_t(5)); ++frame)
	{
		std::cout << "  Frame " << frame << ": (" 
		          << history[0][frame].X() << ", " 
		          << history[0][frame].Y() << ", " 
		          << history[0][frame].Z() << ")" << std::endl;
	}

	// Prepare visualization data
	std::vector<std::string> colors(numParticles, color);
	Box3D bbox = sim.getBoundingBox();
	Real charSize = std::max({bbox.Width(), bbox.Height(), bbox.Depth()});
	std::vector<Real> radii(numParticles, 0.02 * charSize * particleSizeScale);

	// Convert history to Pnt3Cart format
	std::vector<std::vector<Pnt3Cart>> positions(numParticles);
	for (int i = 0; i < numParticles; ++i)
	{
		for (const auto& pt : history[i])
		{
			positions[i].push_back(Pnt3Cart(pt.X(), pt.Y(), pt.Z()));
		}
	}

	// Save and visualize
	std::string fullPath = GetResultFilesPath() + fileName;
	Real boxX = 1.5 * bbox.Width();
	Real boxY = 1.5 * bbox.Height();
	Real boxZ = 1.5 * bbox.Depth();
	Real dt = charSize * 0.0001 * saveEveryN;

	auto result = Serializer::SaveParticleSimulation3D(
		fullPath, numParticles, boxX, boxY, boxZ,
		positions, colors, radii, dt
	);

	if (result.success)
	{
		std::cout << "Saved simulation to: " << fullPath << std::endl;
		std::cout << "Recorded " << positions[0].size() << " frames" << std::endl;
		
		std::cout << "\nLaunching visualizer..." << std::endl;
		std::string vizPath = "tools\\visualizers\\win\\WPF\\particle_3d_visualizer\\MML_ParticleVisualizer3D.exe";
		std::string command = "start \"\" \"" + vizPath + "\" \"" + fullPath + "\"";
		int ret = std::system(command.c_str());
		if (ret != 0) {
			std::cout << "Note: Visualizer may have failed to launch (code " << ret << ")" << std::endl;
			std::cout << "Try running manually:" << std::endl;
			std::cout << "  .\\" << vizPath << " \"" << fullPath << "\"" << std::endl;
		}
	}
	else
	{
		std::cout << "Error saving simulation: " << result.message << std::endl;
	}
}

void Demo1_SphereGeometry(Real radius, int numParticles)
{
	// Create sphere geometry centered at origin
	auto sphere = std::make_unique<Sphere3D>(radius);
	GenericChargeSimulator sim(std::move(sphere));
	
	std::string fileName = "charge_distribution_sphere_R" + std::to_string(static_cast<int>(radius)) + 
	                       "_N" + std::to_string(numParticles) + ".mml";
	runChargeSimulationDemo(sim, "Sphere (R=" + std::to_string(static_cast<int>(radius)) + ")", 
	                        fileName, "Blue", numParticles, 0.2);  // 20% size
}

void Demo1_CubeGeometry(Real sideLen, int numParticles)
{
	// Create cube geometry centered at origin
	auto cube = std::make_unique<Cube3D>(sideLen);
	GenericChargeSimulator sim(std::move(cube));
	
	std::string fileName = "charge_distribution_cube_S" + std::to_string(static_cast<int>(sideLen)) + 
	                       "_N" + std::to_string(numParticles) + ".mml";
	runChargeSimulationDemo(sim, "Cube (side=" + std::to_string(static_cast<int>(sideLen)) + ")",
	                        fileName, "Red", numParticles, 0.3);  // 30% size
}

void Demo1_TubeGeometry(Real radius, Real length, int numParticles)
{
	// Create cylinder geometry centered at origin
	// Note: Cylinder3D constructor takes (radius, height, bottom_center)
	// For centered cylinder, we set bottom at z = -length/2
	auto cylinder = std::make_unique<Cylinder3D>(radius, length, Pnt3Cart(0, 0, -length/2));
	GenericChargeSimulator sim(std::move(cylinder));
	
	std::string fileName = "charge_distribution_cylinder_R" + std::to_string(static_cast<int>(radius)) + 
	                       "_L" + std::to_string(static_cast<int>(length)) +
	                       "_N" + std::to_string(numParticles) + ".mml";
	runChargeSimulationDemo(sim, "Cylinder (R=" + std::to_string(static_cast<int>(radius)) + 
	                        ", L=" + std::to_string(static_cast<int>(length)) + ")",
	                        fileName, "Green", numParticles, 0.2);  // 20% size
}

void Chapter19_electric_charge_distribution()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****          EXAMPLE 17 - electric charge distributions           ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "\nUsing GenericChargeSimulator with MML's IBody geometry interface!" << std::endl;
	std::cout << "Same physics code works for Sphere3D, Cube3D, Cylinder3D, and any IBody!" << std::endl;

	// Demo 1: Charges on a sphere surface (Thomson problem-like)
	// Classic Thomson problem: find minimum energy configuration of N equal charges on a sphere
	// For known solutions:
	//   N=2:  Antipodal points
	//   N=3:  Equilateral triangle
	//   N=4:  Tetrahedron
	//   N=6:  Octahedron
	//   N=8:  Square antiprism
	//   N=12: Icosahedron
	//   N=20: Dodecahedral-like

	Demo1_SphereGeometry(1000.0, 2000);  // Radius 1000, 2000 particles
	
	Demo1_CubeGeometry(1000.0, 2000);    // Side 1000, 2000 particles
	
	Demo1_TubeGeometry(300.0, 1000.0, 2000);  // Radius 300, Length 1000, 2000 particles
}

