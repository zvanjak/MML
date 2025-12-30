#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "mml/MMLBase.h"

#include "mml/base/BaseUtils.h"
#include "mml/base/Vector.h"

#include "mml/tools/ConsolePrinter.h"
#include "mml/tools/Visualizer.h"
#include "mml/tools/Serializer.h"

#include "mml/algorithms/ODESystemSolver.h"
#include "mml/algorithms/ODESystemStepCalculators.h"
#include "mml/algorithms/ODESystemSteppers.h"

#include "mpl/Mechanics/Projectiles2D.h"
#endif

#include <string>
#include <vector>

using namespace MML;
using namespace MPL;

void Projectile_trajectory(IProjectile2DODE &projSys, Real angle, Real initVel, Real initHeight, std::string title)
{
	ProjectileMotionSolver2D projSolver(projSys);

	Real t2 = 1.1 * projSolver.getTimeOfFlightVacuum(angle, initHeight, initVel); 
 
	auto results = projSolver.solveEuler(angle, initHeight, initVel, t2, 0.5);

	// Visualizing the solution
	LinearInterpRealFunc	sol_X_of_t = results.getXOfT();
	LinearInterpRealFunc 	sol_Y_of_t = results.getYOfT();
	LinearInterpRealFunc	sol_Vx_of_t = results.getVxOfT();
	LinearInterpRealFunc	sol_Vy_of_t = results.getVyOfT();

	Visualizer::VisualizeRealFunction(sol_Y_of_t, title + " - height in time",
																		0.0, t2, 500, title + "_height_in_time.txt");
	Visualizer::VisualizeRealFunction(sol_Vx_of_t, title + " - vx in time",
																		0.0, t2, 500, title + "_vx_in_time.txt");
	Visualizer::VisualizeRealFunction(sol_Vy_of_t, title + " - vy in time",
																		0.0, t2, 500, title + "_vy_in_time.txt");

	LinearInterpRealFunc sol_Y_of_X = results.getYOfX();	// forming y(x) function

	Real x2 = 1.1 * results._range;
	Visualizer::VisualizeRealFunction(sol_Y_of_X, title + " - y(x)", 0, x2, 500, title + "_y(x).txt");

	Visualizer::VisualizeODESysSolAsMultiFunc(projSolver.getLastODESolution(), title + " - all variables as func. of t",
																						std::vector<std::string>{"x", "y", "vx", "vy"}, title + "_all_variables.txt");
}

void Projectile_trajectory_multiple_angles(IProjectile2DODE& projSys, Vector<Real> angles, Real initVel, Real initHeight, std::string title)
{
	ProjectileMotionSolver2D projSolver(projSys);

	Vector<ProjectileTrajectory2D> results = projSolver.solveForAnglesEuler(angles, initHeight, initVel, 0.5);

	std::vector<LinearInterpRealFunc> solutions;
	std::vector<std::string> strAngles;
	Real maxRange = 0.0;								// to find the maximum range for visualization
	for (auto& result : results)
	{
		if (result._range > maxRange)
			maxRange = result._range;
		solutions.push_back(result.getYOfX());
		strAngles.push_back("Angle - " + std::to_string(Utils::RadToDeg(result._angle)));
	}

	// now visualize all solutions
	Visualizer::VisualizeMultiRealFunction(solutions,
																				title + " - different angles", strAngles,
																				0, 1.1 * maxRange, 500, title + "_different_angles.txt");
}

void Compare_solutions_for_set_of_angles(IProjectile2DODE& sys1, std::string label1, IProjectile2DODE& sys2, std::string label2,
																				 Real initVel, Real initHeight, Real dt, Vector<Real> angles)
{
	ProjectileMotionSolver2D projSolver1(sys1);
	ProjectileMotionSolver2D projSolver2(sys2);

	Vector<ProjectileTrajectory2D> results1 = projSolver1.solveForAnglesEuler(angles, initHeight, initVel, dt);
	Vector<ProjectileTrajectory2D> results2 = projSolver2.solveForAnglesEuler(angles, initHeight, initVel, dt);
	
	std::vector<LinearInterpRealFunc> solutions;
	std::vector<std::string> strAngles;
	Real maxRange = 0.0;								
	for (int i = 0; i < results1.size(); i++)
	{
		if (results1[i]._range > maxRange)
			maxRange = results1[i]._range;
		solutions.push_back(results1[i].getYOfX());
		strAngles.push_back(label1 + " - " + std::to_string(Utils::RadToDeg(results1[i]._angle)));
	}
	for(int i = 0; i < results2.size(); i++)
	{
		if (results2[i]._range > maxRange)
			maxRange = results2[i]._range;
		solutions.push_back(results2[i].getYOfX());
		strAngles.push_back(label2 + " - " + std::to_string(Utils::RadToDeg(results2[i]._angle)) + " (" + label2 + ")");
	}

	// now visualize all solutions
	Visualizer::VisualizeMultiRealFunction(solutions,
																				"Comparing " + label1 + " and " + label2 + " solutions", strAngles,
																				0, 1.1 * maxRange, 500, label1 + "_multi.txt");

}
