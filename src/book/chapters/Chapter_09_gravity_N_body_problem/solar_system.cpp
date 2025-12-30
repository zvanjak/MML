#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/tools/Serializer.h"
#include "mml/tools/Visualizer.h"

#include "mpl/Base/SolarSystem.h"
#include "mpl/Gravity/NBodySimulator.h"
#endif

using namespace MML;
using namespace MPL;

// demo of Solar system simulation
void Demo_Solar_system()
{
	NBodyGravitySimConfig config = NBodyGravityConfigGenerator::Config2_Solar_system();

	NBodyGravitySimulator solver(config);

	Real t1 = 0.0, t2 = 2;
	const int  steps = 501;
	const Real dt = (t2 - t1) / steps;

	// solving with Euler method
	NBodyGravitySimulationResults result = solver.SolveEuler(dt, steps);

	result.VisualizeAsParamCurve("solar_system", Vector<int>{ 1, 2, 3, 4, 5, 6, 7 });

	// visualize as particle simulation
	result.VisualizeAsParticleSimulation("solar_system_particles", Vector<int>{ 1, 2, 3, 4, 5 }, dt);
}

// demo of star clusters collision
void Demo_star_clusters_collision()
{
	// TODO: Implementation pending
}
