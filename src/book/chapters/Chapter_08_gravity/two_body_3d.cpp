#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/VectorTypes.h"

#include "mml/tools/Serializer.h"
#include "mml/tools/Visualizer.h"
#include "mml/tools/ConsolePrinter.h"

#include "mpl/Gravity/TwoBodySimulator.h"
#endif

#include <iostream>

using namespace MML;
using namespace MPL;

// Two bodies in 3D
// Check conservation of momentum, energy and angular momentum
void Demo_TwoMasses3D()
{
	//TwoBodyGravitySimConfig config = TwoBodyGravityConfigGenerator::Config1_same_bodies_elliptic_CM_static();
	TwoBodyGravitySimConfig config = TwoBodyGravityConfigGenerator::Config2_same_bodies_elliptic_CM_moving();

	TwoBodiesGravitySimulator sim(config);

	auto result = sim.SolveRK5(300, 1e-08, 0.1);

	std::cout << "Results for two bodies simulation in 3D" << std::endl;
	std::cout << "t = " << result._duration << std::endl;
	std::cout << "Mass1 = " << result._config.Mass1() << std::endl;
	std::cout << "Mass2 = " << result._config.Mass2() << std::endl;
	std::cout << "G = " << result._config.G() << std::endl;
	std::cout << "Number of steps = " << result._vecTimes.size() << std::endl;

	std::cout << "\nPosition and velocity values for both bodies" << std::endl;
	int fmt1W = 10, fmt1D = 5;
	std::vector<ColumnFormat> vecNames{ ColumnFormat("t", 11, 2, 'F'),
																 ColumnFormat("x1", fmt1W, fmt1D, 'F'), ColumnFormat("y1", fmt1W, fmt1D, 'F'), ColumnFormat("z1", fmt1W, fmt1D, 'F'),
																 ColumnFormat("x2", fmt1W, fmt1D, 'F'), ColumnFormat("y2", fmt1W, fmt1D, 'F'), ColumnFormat("z2", fmt1W, fmt1D, 'F'),
																 ColumnFormat("v1_x", fmt1W, fmt1D, 'F'), ColumnFormat("v1_y", fmt1W, fmt1D, 'F'), ColumnFormat("v1_z", fmt1W, fmt1D, 'F'),
																 ColumnFormat("v2_x", fmt1W, fmt1D, 'F'), ColumnFormat("v2_y", fmt1W, fmt1D, 'F'), ColumnFormat("v2_z", fmt1W, fmt1D, 'F') };

	Vector<Real> t_vals = result.getTimes();
	Vector<Real> x1_vals = result.getPos1X();
	Vector<Real> y1_vals = result.getPos1Y();
	Vector<Real> z1_vals = result.getPos1Z();
	Vector<Real> x2_vals = result.getPos2X();
	Vector<Real> y2_vals = result.getPos2Y();
	Vector<Real> z2_vals = result.getPos2Z();

	Vector<Real> v1_x_vals = result.getV1X();
	Vector<Real> v1_y_vals = result.getV1Y();
	Vector<Real> v1_z_vals = result.getV1Z();
	Vector<Real> v2_x_vals = result.getV2X();
	Vector<Real> v2_y_vals = result.getV2Y();
	Vector<Real> v2_z_vals = result.getV2Z();

	std::vector<Vector<Real>*> vecVals{ &t_vals, &x1_vals, &y1_vals, &z1_vals, &x2_vals, &y2_vals, &z2_vals,
																			&v1_x_vals,&v1_y_vals,&v1_z_vals,&v2_x_vals,&v2_y_vals,&v2_z_vals };

	VectorTablePrinter vvp(vecNames, vecVals);
	vvp.print();

	std::cout << "\nCenter of mass position and velocity values" << std::endl;
	int fmt2W = 10, fmt2D = 5;
	std::vector<ColumnFormat> vecNames2{ ColumnFormat("t", 11, 2, 'F'),
																	ColumnFormat("cm_x", fmt2W, fmt2D, 'F'), ColumnFormat("cm_y", fmt2W, fmt2D, 'F'), ColumnFormat("cm_z", fmt2W, fmt2D, 'F'),
																	ColumnFormat("cm_vx", fmt2W, fmt2D, 'F'), ColumnFormat("cm_vy", fmt2W, fmt2D, 'F'), ColumnFormat("cm_vz", fmt2W, fmt2D, 'F') };

	Vector<Real> cm_x_vals = result.getCMPosX();
	Vector<Real> cm_y_vals = result.getCMPosY();
	Vector<Real> cm_z_vals = result.getCMPosZ();

	Vector<Real> cm_vx_vals = result.getCMVX();
	Vector<Real> cm_vy_vals = result.getCMVY();
	Vector<Real> cm_vz_vals = result.getCMVZ();

	std::vector<Vector<Real>*> vecVals2{ &t_vals, &cm_x_vals, &cm_y_vals, &cm_z_vals,
																			 &cm_vx_vals, &cm_vy_vals, &cm_vz_vals };

	VectorTablePrinter vvp2(vecNames2, vecVals2);
	vvp2.print();

	std::cout << "\nTotal momentum and energy values" << std::endl;
	int fmt3W = 12, fmt3D = 2;
	std::vector<ColumnFormat> vecNames3{ ColumnFormat("t", 11, 2, 'F'),
																	ColumnFormat("p_x", fmt3W, fmt3D, 'F'), ColumnFormat("p_y", fmt3W, fmt3D, 'F'), ColumnFormat("p_z", fmt3W, fmt3D, 'F'),
																	ColumnFormat("ke", fmt3W, fmt3D, 'F'), ColumnFormat("pe", fmt3W, fmt3D, 'F'), ColumnFormat("te", fmt3W, fmt3D, 'F') };

	Vector<Real> p_x_vals = result.getTotalMomentumX();
	Vector<Real> p_y_vals = result.getTotalMomentumY();
	Vector<Real> p_z_vals = result.getTotalMomentumZ();

	Vector<Real> ke_vals = result.getTotalKineticEnergy();
	Vector<Real> pe_vals = result.getTotalPotentialEnergy();
	Vector<Real> te_vals = result.getTotalEnergy();

	std::vector<Vector<Real>*> vecVals3{ &t_vals, &p_x_vals, &p_y_vals, &p_z_vals,
																			 &ke_vals, &pe_vals, &te_vals };

	VectorTablePrinter vvp3(vecNames3, vecVals3);
	vvp3.print();

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
	vvp4.print();

	result.VisualizeSolution("gravity_example_3d");
}
