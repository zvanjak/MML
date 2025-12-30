#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Derivation.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/ODESystemSteppers.h"

#include "tools/ConsolePrinter.h"
#include "tools/Visualizer.h"
#include "tools/Serializer.h"

#include "mpl/Mechanics/Pendulums.h"
#endif

using namespace MML;
using namespace MPL;

void Demo_DoublePendulum()
{
	Real	l1 = 1.0, l2 = 1.0;
	Real  m1 = 0.5, m2 = 1.0;
	DoublePendulumODE		odeSysDoublePend = DoublePendulumODE(m1, m2, l1, l2);

	Real	t1 = 0.0, t2 = 50.0;
	int   expectNumSteps = 10000;
	Real	minSaveInterval = (t2 - t1) / expectNumSteps;
	Real	initAngle1 = 0.5;
	Real  initAngle2 = 0.101;
	Vector<Real>	initCond{ initAngle1, 0.0, initAngle2, 0.0 };

	ODESystemSolver<RK5_CashKarp_Stepper> solver(odeSysDoublePend);
	ODESystemSolution sol = solver.integrate(initCond, t1, t2, minSaveInterval, 1e-08, 0.01);

	Vector<Real> t_vals = sol.getTValues();
	Vector<Real> theta1_vals = sol.getXValues(0);
	Vector<Real> theta2_vals = sol.getXValues(2);

	std::cout << "\n\n**** Solving double pendulum  ****\n";
	std::vector<ColDesc>				vecNames{ ColDesc("t", 11, 2, 'F'),
																				ColDesc("theta 1", 15, 8, 'F'),
																				ColDesc("theta 2", 15, 8, 'F'), };
	std::vector<Vector<Real>*>	vecVals{ &t_vals, &theta1_vals, &theta2_vals };

	VectorTablePrinter	vvp(vecNames, vecVals);

	//vvpprint();

	// getting solutions as polynomials
	PolynomInterpRealFunc 	solAdaptPolyInterp0 = sol.getSolAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solAdaptPolyInterp1 = sol.getSolAsPolyInterp(2, 3);

	Visualizer::VisualizeRealFunction(solAdaptPolyInterp0,
		"Double pendulum - theta1 in time",
		0.0, 10.0, 200, "double_pendulum_angle1.txt");

	// shown together
	Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{ &solAdaptPolyInterp0, &solAdaptPolyInterp1 },
		"Double pendulum - both angles", { "Angle", "Ang.vel" },
		0.0, 10.0, 200,
		"double_pendulum_multi_real_func.txt");

	// form parametric curve 2d from solution
	Matrix<Real> curve_points(t_vals.size(), 2);
	for (int i = 0; i < t_vals.size(); i++)
	{
		curve_points(i, 0) = theta1_vals[i];
		curve_points(i, 1) = theta2_vals[i];
	}

	SplineInterpParametricCurve<2> phaseSpaceTrajectory(0.0, 1.0, curve_points);

	Visualizer::VisualizeParamCurve2D(phaseSpaceTrajectory, "Double pendulum - phase space trajectory",
		0.0, 1.0, t_vals.size(), "double_pendulum_phase_space.txt");
}


void Demo_SphericalPendulum()
{
	Real	l = 1.0;
	//SphericalPendulumODE		odeSysSpherPend = SphericalPendulumODE(l);
	SphericalPendulumHamiltonODE		odeSysSpherPend = SphericalPendulumHamiltonODE(l);

	Real	t1 = 0.0, t2 = 20.0;
	int   expectNumSteps = 2000;
	Real	minSaveInterval = (t2 - t1) / expectNumSteps;

	Real	initAngleTheta = Constants::PI / 2;
	Real  initAnglePhi = 0.0;
	Vector<Real>	initCond{ initAngleTheta, initAnglePhi, 5.10, -3.55 };

	ODESystemSolver<RK5_CashKarp_Stepper> odeSolver(odeSysSpherPend);
	ODESystemSolution odeSol = odeSolver.integrate(initCond, t1, t2, minSaveInterval, 1e-08, 0.01);

	Vector<Real> t_vals = odeSol.getTValues();
	Vector<Real> theta_vals = odeSol.getXValues(0);
	Vector<Real> phi_vals = odeSol.getXValues(1);
	Vector<Real> theta_dot_vals = odeSol.getXValues(2);
	Vector<Real> phi_dot_vals = odeSol.getXValues(3);

	std::cout << "\n\n****  Runge-Kutta 4th order - fixed stepsize  **********  Runge-Kutta 4th order - adaptive stepper  ****\n";
	std::vector<ColDesc>				vecNames{ ColDesc("t", 11, 2, 'F'), ColDesc("theta", 15, 8, 'F'), ColDesc("theta dot", 15, 8, 'F'),
																				ColDesc("phi", 12, 8, 'F'), ColDesc("phi dot", 12, 8, 'F') };
	std::vector<Vector<Real>*>	vecVals{ &t_vals, &theta_vals, &theta_dot_vals, &phi_vals, &phi_dot_vals };
	VectorTablePrinter	vvp(vecNames, vecVals);

	//vvpprint();

	// getting solutions as polynomials
	PolynomInterpRealFunc 	solAdaptPolyInterp0 = odeSol.getSolAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solAdaptPolyInterp1 = odeSol.getSolAsPolyInterp(1, 3);

	Visualizer::VisualizeRealFunction(solAdaptPolyInterp0, "Spherical pendulum - theta in time",
		t1, t2, 1000, "spherical_pendulum_theta.txt");

	// shown together
	Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{ &solAdaptPolyInterp0, &solAdaptPolyInterp1 },
		"Spherical pendulum - both angles", { "Theta", "Phi" },
		t1, t2, 1000,
		"spherical_pendulum_multi_real_func.txt");

	// projecting pendulum path on x-y, x-z and y-z plane
	// first, get values for x, y, z coordinates in vectors
	Vector<Real> x_vals(t_vals.size()), y_vals(t_vals.size()), z_vals(t_vals.size());

	for (int i = 0; i < t_vals.size(); i++)
	{
		x_vals[i] = l * sin(theta_vals[i]) * cos(phi_vals[i]);
		y_vals[i] = l * sin(theta_vals[i]) * sin(phi_vals[i]);
		z_vals[i] = l * (1 - cos(theta_vals[i]));
	}

	// then form appropriate matrices with relevant data for visualization of parametric curves
	Matrix<Real> curve_x_y_points(t_vals.size(), 2), curve_x_z_points(t_vals.size(), 2);
	Matrix<Real> curve_y_z_points(t_vals.size(), 2);
	Matrix<Real> curve_xyz_points(t_vals.size(), 3);
	for (int i = 0; i < t_vals.size(); i++)
	{
		curve_xyz_points(i, 0) = 100 * x_vals[i];
		curve_xyz_points(i, 1) = 100 * y_vals[i];
		curve_xyz_points(i, 2) = 100 * z_vals[i];

		curve_x_y_points(i, 0) = x_vals[i];
		curve_x_y_points(i, 1) = y_vals[i];

		curve_x_z_points(i, 0) = x_vals[i];
		curve_x_z_points(i, 1) = z_vals[i];

		curve_y_z_points(i, 0) = y_vals[i];
		curve_y_z_points(i, 1) = z_vals[i];
	}

	SplineInterpParametricCurve<3> xyzTrajectory(0.0, 1.0, curve_xyz_points);
	Visualizer::VisualizeParamCurve3D(xyzTrajectory, "Spherical pendulum : x-y-z trajectory",
		0.0, 1.0, t_vals.size(), "spherical_pendulum_xyz_trajectory.txt");

	SplineInterpParametricCurve<2> xyTrajectory(0.0, 1.0, curve_x_y_points);
	Visualizer::VisualizeParamCurve2D(xyTrajectory, "Spherical pendulum : x-y trajectory",
		0.0, 1.0, t_vals.size(), "spherical_pendulum_xy_trajectory.txt");

	SplineInterpParametricCurve<2> xzTrajectory(0.0, 1.0, curve_x_z_points);
	Visualizer::VisualizeParamCurve2D(xzTrajectory, "Spherical pendulum : x-z trajectory",
		0.0, 1.0, t_vals.size(), "spherical_pendulum_xz_trajectory.txt");

	SplineInterpParametricCurve<2> yzTrajectory(0.0, 1.0, curve_y_z_points);
	Visualizer::VisualizeParamCurve2D(yzTrajectory, "Spherical pendulum : y-z trajectory",
		0.0, 1.0, t_vals.size(), "spherical_pendulum_yz_trajectory.txt");
}

void Example7_Lagrangian_mechanics()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                   EXAMPLE 6 - Lagrangian mechanics            ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	//Demo_DoublePendulum();
	Demo_SphericalPendulum();
}