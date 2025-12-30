#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/VectorTypes.h"
#include "base/Function.h"

#include "core/Curves.h"
#include "core/Derivation.h"
#include "core/Integration/PathIntegration.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/ODESystemSteppers.h"

#include "tools/Serializer.h"
#include "tools/Visualizer.h"
#include "tools/ConsolePrinter.h"

#include "mpl/Gravity/GravityBase.h"
#include "mpl/Gravity/TwoBodySimulator.h"
#include "mpl/Gravity/TwoBodySimulator2D.h"
#endif

using namespace MML;
using namespace MPL;

// calculate, simulate and visualize motion for two bodies in 2D
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
	std::vector<ColDesc> vecNames{ ColDesc("t", 11, 2, 'F'),
																 ColDesc("x1", fmt1W, fmt1D, 'F'), ColDesc("y1", fmt1W, fmt1D, 'F'), 
																 ColDesc("x2", fmt1W, fmt1D, 'F'), ColDesc("y2", fmt1W, fmt1D, 'F'), 
																 ColDesc("v1_x", fmt1W, fmt1D, 'F'), ColDesc("v1_y", fmt1W, fmt1D, 'F'), 
																 ColDesc("v2_x", fmt1W, fmt1D, 'F'), ColDesc("v2_y", fmt1W, fmt1D, 'F') };

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
	vvpprint();

	std::cout << "\nCenter of mass position and velocity values" << std::endl;
	int fmt2W = 10, fmt2D = 5;
	std::vector<ColDesc> vecNames2{ ColDesc("t", 11, 2, 'F'),
																	ColDesc("cm_x", fmt2W, fmt2D, 'F'), ColDesc("cm_y", fmt2W, fmt2D, 'F'), 
																	ColDesc("cm_vx", fmt2W, fmt2D, 'F'), ColDesc("cm_vy", fmt2W, fmt2D, 'F') };

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
	std::vector<ColDesc> vecNames3{ ColDesc("t", 11, 2, 'F'),
																	ColDesc("p_x", fmt3W, fmt3D, 'F'), ColDesc("p_y", fmt3W, fmt3D, 'F'),
																	ColDesc("ke", fmt3W, fmt3D, 'F'), ColDesc("pe", fmt3W, fmt3D, 'F'), ColDesc("te", fmt3W, fmt3D, 'F') };

	Vector<Real> p_x_vals = result.getTotalMomentumX();
	Vector<Real> p_y_vals = result.getTotalMomentumY();

	Vector<Real> ke_vals = result.getTotalKineticEnergy();
	Vector<Real> pe_vals = result.getTotalPotentialEnergy();
	Vector<Real> te_vals = result.getTotalEnergy();

	std::vector<Vector<Real>*> vecVals3{ &t_vals, &p_x_vals, &p_y_vals, 
																			 &ke_vals, &pe_vals, &te_vals };

	VectorTablePrinter vvp3(vecNames3, vecVals3);
	vvp3print();



}

// checking Kepler laws
void Demo_Check_Kepler_laws()
{

}

// two bodies in 3D
// check conservation of momentum, energy and angular momentum
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
	std::vector<ColDesc> vecNames{ ColDesc("t", 11, 2, 'F'),
																 ColDesc("x1", fmt1W, fmt1D, 'F'), ColDesc("y1", fmt1W, fmt1D, 'F'), ColDesc("z1", fmt1W, fmt1D, 'F'),
																 ColDesc("x2", fmt1W, fmt1D, 'F'), ColDesc("y2", fmt1W, fmt1D, 'F'), ColDesc("z2", fmt1W, fmt1D, 'F'),
																 ColDesc("v1_x", fmt1W, fmt1D, 'F'), ColDesc("v1_y", fmt1W, fmt1D, 'F'), ColDesc("v1_z", fmt1W, fmt1D, 'F'),
																 ColDesc("v2_x", fmt1W, fmt1D, 'F'), ColDesc("v2_y", fmt1W, fmt1D, 'F'), ColDesc("v2_z", fmt1W, fmt1D, 'F') };

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
	vvpprint();

	std::cout << "\nCenter of mass position and velocity values" << std::endl;
	int fmt2W = 10, fmt2D = 5;
	std::vector<ColDesc> vecNames2{ ColDesc("t", 11, 2, 'F'),
																	ColDesc("cm_x", fmt2W, fmt2D, 'F'), ColDesc("cm_y", fmt2W, fmt2D, 'F'), ColDesc("cm_z", fmt2W, fmt2D, 'F'),
																	ColDesc("cm_vx", fmt2W, fmt2D, 'F'), ColDesc("cm_vy", fmt2W, fmt2D, 'F'), ColDesc("cm_vz", fmt2W, fmt2D, 'F') };

	Vector<Real> cm_x_vals = result.getCMPosX();
	Vector<Real> cm_y_vals = result.getCMPosY();
	Vector<Real> cm_z_vals = result.getCMPosZ();

	Vector<Real> cm_vx_vals = result.getCMVX();
	Vector<Real> cm_vy_vals = result.getCMVY();
	Vector<Real> cm_vz_vals = result.getCMVZ();

	std::vector<Vector<Real>*> vecVals2{ &t_vals, &cm_x_vals, &cm_y_vals, &cm_z_vals,
																			 &cm_vx_vals, &cm_vy_vals, &cm_vz_vals };

	VectorTablePrinter vvp2(vecNames2, vecVals2);
	vvp2print();

	std::cout << "\nTotal momentum and energy values" << std::endl;
	int fmt3W = 12, fmt3D = 2;
	std::vector<ColDesc> vecNames3{ ColDesc("t", 11, 2, 'F'),
																	ColDesc("p_x", fmt3W, fmt3D, 'F'), ColDesc("p_y", fmt3W, fmt3D, 'F'), ColDesc("p_z", fmt3W, fmt3D, 'F'),
																	ColDesc("ke", fmt3W, fmt3D, 'F'), ColDesc("pe", fmt3W, fmt3D, 'F'), ColDesc("te", fmt3W, fmt3D, 'F') };

	Vector<Real> p_x_vals = result.getTotalMomentumX();
	Vector<Real> p_y_vals = result.getTotalMomentumY();
	Vector<Real> p_z_vals = result.getTotalMomentumZ();

	Vector<Real> ke_vals = result.getTotalKineticEnergy();
	Vector<Real> pe_vals = result.getTotalPotentialEnergy();
	Vector<Real> te_vals = result.getTotalEnergy();

	std::vector<Vector<Real>*> vecVals3{ &t_vals, &p_x_vals, &p_y_vals, &p_z_vals,
																			 &ke_vals, &pe_vals, &te_vals };

	VectorTablePrinter vvp3(vecNames3, vecVals3);
	vvp3print();

	std::cout << "\nAngular momentum values" << std::endl;
	int fmt4W = 12, fmt4D = 2;
	std::vector<ColDesc> vecNames4{ ColDesc("t", 11, 2, 'F'),
																	ColDesc("L_x", fmt4W, fmt4D, 'F'), ColDesc("L_y", fmt4W, fmt4D, 'F'), ColDesc("L_z", fmt4W, fmt4D, 'F'),
																	ColDesc("L_x_orig", fmt4W, fmt4D, 'F'), ColDesc("L_y_orig", fmt4W, fmt4D, 'F'), ColDesc("L_z_orig", fmt4W, fmt4D, 'F') };

	Vector<Real> L_x_vals = result.getAngularMomentumCMX();
	Vector<Real> L_y_vals = result.getAngularMomentumCMY();
	Vector<Real> L_z_vals = result.getAngularMomentumCMZ();

	Vector<Real> L_x_orig_vals = result.getAngularMomentumX(Vec3Cart(0, 0, 0));
	Vector<Real> L_y_orig_vals = result.getAngularMomentumY(Vec3Cart(0, 0, 0));
	Vector<Real> L_z_orig_vals = result.getAngularMomentumZ(Vec3Cart(0, 0, 0));

	std::vector<Vector<Real>*> vecVals4{ &t_vals, &L_x_vals, &L_y_vals, &L_z_vals,
																			 &L_x_orig_vals, &L_y_orig_vals, &L_z_orig_vals };
	VectorTablePrinter vvp4(vecNames4, vecVals4);
	vvp4print();

	result.VisualizeSolution("gravity_example_3d");
}

// verify vector and scalar field operations for gravity fields
void Demo_Field_operations()
{
}

// for two given points in space, calculates the potential at those points, 
// and then calculates the path integral of the force field along a line and a spline curve bdsetween those two points
void Verify_path_integrals()
{
	Real G = 1.0;
	Real mass = 100.0;

	GravityPotentialField		potentialField(G, mass, Vec3Cart(0, 0, 0));
	GravityForceField				forceField(G, mass, Vec3Cart(0, 0, 0));

	// our two points in space
	Vec3Cart point1(10.0, 10.0, -5.0);
	Vec3Cart point2(-5.0, 5.0, 8.0);

	Real potential1 = potentialField(point1);
	Real potential2 = potentialField(point2);

	std::cout << "Potential at point 1: " << potential1 << std::endl;
	std::cout << "Potential at point 2: " << potential2 << std::endl;
	std::cout << "Difference in potential: " << potential1 - potential2 << std::endl;

	// first, create a line between those two points
	Curves::LineCurve line(0.0, Pnt3Cart(point1.X(), point1.Y(), point1.Z()), 
												 1.0, Pnt3Cart(point2.X(), point2.Y(), point2.Z()));

	Real integral = PathIntegration::LineIntegral(forceField, line, 0.0, 1.0);

	std::cout << "Path integral between point 1 and point 2: " << integral << std::endl;

	// second, we'll create a spline parametric curve between those two points
	// with some points added in between to make path "curvaceus"
	Matrix<Real> curve_points1{ 5, 3,
														 { point1.X(), point1.Y(), point1.Z(),
															20.0, 5.0, -10.0,
															100.0, 150.5, -30.0,		// sending it far away
															5.0, 10.5, 1.0,
															point2.X(), point2.Y(), point2.Z() }
	};
	SplineInterpParametricCurve<3> curve(0.0, 1.0, curve_points1);

	Real integral2 = PathIntegration::LineIntegral(forceField, curve, 0.0, 1.0);

	std::cout << "Path integral along spline curve between point 1 and point 2: " << integral2 << std::endl;
}

// simulate gravitational slingshot of Voyger 1 by Jupiter
void Demo_Voyager_Jupiter_Slingshot()
{
	// to be implemented
}

void Example8_Gravity()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                       EXAMPLE 8 - Gravity                     ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Demo_TwoMasses2D();
	//Demo_Check_Kepler_laws();
	//Demo_TwoMasses3D();
	//Demo_Field_operations();
	//Verify_path_integrals();
	//Demo_Voyager_Jupiter_Slingshot();
}