#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"


#include "mpl/Base/SolarSystem.h"
#include "mpl/Gravity/NBodySimulator.h"

#include "mml/tools/Serializer.h"
#include "mml/tools/Visualizer.h"
#include "mml/tools/ConsolePrinter.h"

#endif

#include <iostream>

using namespace MML;
using namespace MPL;

// demo with 5 masses, with different initial positions and velocities
void Demo_FiveMasses()
{
	NBodyGravitySimConfig config = NBodyGravityConfigGenerator::Config1_five_bodies();

	NBodyGravitySimulator solver(config);

	Real t1 = 0.0, t2 = 200.0;

	// solving with Euler method
	const int  steps = 2001;
	const Real dt = (t2 - t1) / steps;
	NBodyGravitySimulationResults result = solver.SolveEuler(dt, steps);

	// VISUALIZING TRAJECTORIES AS PARAMETRIC CURVES
	result.VisualizeAsParamCurve("five_bodies", Vector<int>{ 0, 1, 2, 3, 4 });

	//NBodyGravitySimulationResults result = solver.SolveRK5(t2, 1e-08, 0.1, 0.1);

	// PRINTING RESULTS TO CONSOLE
	std::cout << "\nCenter of mass position and velocity values" << std::endl;
	int fmt2W = 12, fmt2D = 7;
	std::vector<ColumnFormat> vecNames2{ ColumnFormat("t", 11, 2, 'F'),
																	ColumnFormat("cm_x", fmt2W, fmt2D, 'F'), ColumnFormat("cm_y", fmt2W, fmt2D, 'F'), ColumnFormat("cm_z", fmt2W, fmt2D, 'F'),
																	ColumnFormat("cm_vx", fmt2W, fmt2D, 'F'), ColumnFormat("cm_vy", fmt2W, fmt2D, 'F'), ColumnFormat("cm_vz", fmt2W, fmt2D, 'F') };
	
	Vector<Real> t_vals = result.getTimes();
	Vector<Real> cm_x_vals = result.getCMPosX();
	Vector<Real> cm_y_vals = result.getCMPosY();
	Vector<Real> cm_z_vals = result.getCMPosZ();

	Vector<Real> cm_vx_vals = result.getCMVX();
	Vector<Real> cm_vy_vals = result.getCMVY();
	Vector<Real> cm_vz_vals = result.getCMVZ();

	std::vector<Vector<Real>*> vecVals2{ &t_vals, &cm_x_vals, &cm_y_vals, &cm_z_vals,
																			 &cm_vx_vals, &cm_vy_vals, &cm_vz_vals };

	VectorTablePrinter vvp2(vecNames2, vecVals2);
	//vvp2print();

	std::cout << "\nTotal momentum and energy values" << std::endl;
	int fmt3W = 12, fmt3D = 2;
	std::vector<ColumnFormat> vecNames3{ ColumnFormat("t", 11, 2, 'F'),
																	ColumnFormat("p_x", fmt3W, fmt3D, 'F'), ColumnFormat("p_y", fmt3W, fmt3D, 'F'), ColumnFormat("p_z", fmt3W, fmt3D, 'F'),
																	ColumnFormat("ke", fmt3W, fmt3D, 'F'), ColumnFormat("pe", fmt3W, fmt3D, 'F'), ColumnFormat("te", fmt3W, fmt3D, 'F') };

	Vector<Real> p_x_vals = result.getLinearMomentumX();
	Vector<Real> p_y_vals = result.getLinearMomentumY();
	Vector<Real> p_z_vals = result.getLinearMomentumZ();

	Vector<Real> ke_vals = result.getTotalKineticEnergy();
	Vector<Real> pe_vals = result.getTotalPotentialEnergy();
	Vector<Real> te_vals = result.getTotalEnergy();

	std::vector<Vector<Real>*> vecVals3{ &t_vals, &p_x_vals, &p_y_vals, &p_z_vals,
																			 &ke_vals, &pe_vals, &te_vals };

	VectorTablePrinter vvp3(vecNames3, vecVals3);
	//vvp3print();

	std::cout << "\nAngular momentum values" << std::endl;
	int fmt4W = 12, fmt4D = 2;
	std::vector<ColumnFormat> vecNames4{ ColumnFormat("t", 11, 2, 'F'),
																	ColumnFormat("L_x", fmt4W, fmt4D, 'F'), ColumnFormat("L_y", fmt4W, fmt4D, 'F'), ColumnFormat("L_z", fmt4W, fmt4D, 'F'),
																	ColumnFormat("L_x_orig", fmt4W, fmt4D, 'F'), ColumnFormat("L_y_orig", fmt4W, fmt4D, 'F'), ColumnFormat("L_z_orig", fmt4W, fmt4D, 'F') };

	Vector<Real> L_x_vals = result.getAngularMomentumCMX();
	Vector<Real> L_y_vals = result.getAngularMomentumCMY();
	Vector<Real> L_z_vals = result.getAngularMomentumCMZ();

	Vector<Real> L_x_orig_vals = result.getAngularMomentumX(Vec3Cart(0, 0, 0));
	Vector<Real> L_y_orig_vals = result.getAngularMomentumY(Vec3Cart(0, 0, 0));
	Vector<Real> L_z_orig_vals = result.getAngularMomentumZ(Vec3Cart(0, 0, 0));

	std::vector<Vector<Real>*> vecVals4{ &t_vals, &L_x_vals, &L_y_vals, &L_z_vals,
																			 &L_x_orig_vals, &L_y_orig_vals, &L_z_orig_vals };
	VectorTablePrinter vvp4(vecNames4, vecVals4);
	//vvp4print();
}
