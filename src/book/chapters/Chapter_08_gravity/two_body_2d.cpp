#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/VectorTypes.h"

#include "mml/tools/Visualizer.h"
#include "mml/tools/ConsolePrinter.h"

#include "mpl/Gravity/TwoBodySimulator2D.h"
#endif

#include <iostream>
#include <string>

using namespace MML;
using namespace MPL;

// Calculate, simulate and visualize motion for two bodies in 2D
void Demo_TwoMasses2D()
{
	//TwoBodyGravitySimConfig2D config = TwoBodyGravityConfigGenerator2D::Config1_same_bodies_elliptic_CM_static();
	//TwoBodyGravitySimConfig2D config = TwoBodyGravityConfigGenerator2D::Config2_same_bodies_elliptic_CM_moving();
	//TwoBodyGravitySimConfig2D config = TwoBodyGravityConfigGenerator2D::Config3_diff_bodies_elliptic_CM_static();
	//TwoBodyGravitySimConfig2D config = TwoBodyGravityConfigGenerator2D::Config4_diff_bodies_elliptic_CM_moving();
	TwoBodyGravitySimConfig2D config = TwoBodyGravityConfigGenerator2D::Config5_same_bodies_hyperbolic_CM_static();
	//TwoBodyGravitySimConfig2D config = TwoBodyGravityConfigGenerator2D::Config6_same_bodies_hyperbolic_CM_moving();
	
	TwoBodiesGravitySimulator2D sim(config);
	
	//auto result = sim.SolveRK5(100, 1e-08, 0.5);
	//auto result = sim.SolveRK5(250, 1e-08, 0.5);
	auto result = sim.SolveRK5(50000, 1e-08, 0.5);
	//auto result = sim.SolveRK5(1000, 1e-08, 0.5);
	
	result.VisualizeSolution("gravity_example_2d");


	std::cout << "Results for two bodies simulation in 2D" << std::endl;
	std::cout << "t = " << result._duration << std::endl;
	std::cout << "Mass1 = " << result._config.Mass1() << std::endl;
	std::cout << "Mass2 = " << result._config.Mass2() << std::endl;
	std::cout << "G = " << result._config.G() << std::endl;
	std::cout << "Number of steps = " << result.NumSteps() << std::endl;
	
	// calculate orbit type
	std::string trajectoryType = TwoBodyGravityCalculator2D::GetTrajectoryType(config.G(), config.InitState().Body1(), config.InitState().Body2());
	std::cout << "Trajectory type: " << trajectoryType << std::endl;

	std::cout << "\nPosition and velocity values for both bodies" << std::endl;
	int fmt1W = 10, fmt1D = 5;
	std::vector<ColumnFormat> vecNames{ ColumnFormat("t", 11, 2, 'F'),
																 ColumnFormat("x1", fmt1W, fmt1D, 'F'), ColumnFormat("y1", fmt1W, fmt1D, 'F'), 
																 ColumnFormat("x2", fmt1W, fmt1D, 'F'), ColumnFormat("y2", fmt1W, fmt1D, 'F'), 
																 ColumnFormat("v1_x", fmt1W, fmt1D, 'F'), ColumnFormat("v1_y", fmt1W, fmt1D, 'F'), 
																 ColumnFormat("v2_x", fmt1W, fmt1D, 'F'), ColumnFormat("v2_y", fmt1W, fmt1D, 'F') };

	Vector<Real> t_vals = result.getTimes();
	Vector<Real> x1_vals = result.getPos1X();
	Vector<Real> y1_vals = result.getPos1Y();
	Vector<Real> x2_vals = result.getPos2X();
	Vector<Real> y2_vals = result.getPos2Y();

	Vector<Real> v1_x_vals = result.getV1X();
	Vector<Real> v1_y_vals = result.getV1Y();
	Vector<Real> v2_x_vals = result.getV2X();
	Vector<Real> v2_y_vals = result.getV2Y();

	std::vector<Vector<Real>*> vecVals{ &t_vals, &x1_vals, &y1_vals, &x2_vals, &y2_vals,
																			&v1_x_vals,&v1_y_vals, &v2_x_vals,&v2_y_vals };

	VectorTablePrinter vvp(vecNames, vecVals);
	vvp.print();
	
	// TODO: Fix - following code has undefined variables fmt2W, fmt2D
	/*
	std::vector<ColumnFormat> vecNames2{ ColumnFormat("t", 11, 2, 'F'),
																	ColumnFormat("cm_x", fmt2W, fmt2D, 'F'), ColumnFormat("cm_y", fmt2W, fmt2D, 'F'), 
																	ColumnFormat("cm_vx", fmt2W, fmt2D, 'F'), ColumnFormat("cm_vy", fmt2W, fmt2D, 'F') };

	Vector<Real> cm_x_vals = result.getCMPosX();
	Vector<Real> cm_y_vals = result.getCMPosY();

	Vector<Real> cm_vx_vals = result.getCMVX();
	Vector<Real> cm_vy_vals = result.getCMVY();

	std::vector<Vector<Real>*> vecVals2{ &t_vals, &cm_x_vals, &cm_y_vals,
																			 &cm_vx_vals, &cm_vy_vals };

	VectorTablePrinter vvp2(vecNames2, vecVals2);
	vvp2print();

	std::cout << "\nTotal momentum and energy values" << std::endl;
	int fmt3W = 12, fmt3D = 2;
	std::vector<ColumnFormat> vecNames3{ ColumnFormat("t", 11, 2, 'F'),
																	ColumnFormat("p_x", fmt3W, fmt3D, 'F'), ColumnFormat("p_y", fmt3W, fmt3D, 'F'),
																	ColumnFormat("ke", fmt3W, fmt3D, 'F'), ColumnFormat("pe", fmt3W, fmt3D, 'F'), ColumnFormat("te", fmt3W, fmt3D, 'F') };

	Vector<Real> p_x_vals = result.getTotalMomentumX();
	Vector<Real> p_y_vals = result.getTotalMomentumY();

	Vector<Real> ke_vals = result.getTotalKineticEnergy();
	Vector<Real> pe_vals = result.getTotalPotentialEnergy();
	Vector<Real> te_vals = result.getTotalEnergy();

	std::vector<Vector<Real>*> vecVals3{ &t_vals, &p_x_vals, &p_y_vals, 
																			 &ke_vals, &pe_vals, &te_vals };

	VectorTablePrinter vvp3(vecNames3, vecVals3);
	vvp3print();
	*/
}

// Checking Kepler laws
void Demo_Check_Kepler_laws()
{
	// TODO: Implementation pending
}
