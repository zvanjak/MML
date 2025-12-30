#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "mml/MMLBase.h"

#include "mml/base/Vector.h"

#include "mml/tools/Visualizer.h"

#include "mml/algorithms/ODESystemSolver.h"
#include "mml/algorithms/ODESystemSteppers.h"

#include "mpl/Mechanics/Pendulums.h"
#endif

using namespace MML;
using namespace MPL;

void Demo_ForcedPendulum()
{
	ForcedPendulumODE forcedPen(1.0, 0.2, 1);

	Real	t1 = 0.0, t2 = 30.0;
	int   expectNumSteps = 3000;
	Real	minSaveInterval = (t2 - t1) / expectNumSteps;
	Real	initAngle = 0.5;
	Vector<Real>	initCond{ initAngle, 0.0 };

	ODESystemSolver<RK5_CashKarp_Stepper> odeSolver(forcedPen);
	ODESystemSolution sol = odeSolver.integrate(initCond, t1, t2, minSaveInterval, 1e-08, 0.001);

	Vector<Real> t_vals{ sol.getTValues() }, y1_fixed{ sol.getXValues(0) }, y2_fixed{ sol.getXValues(1) };
	
	PolynomInterpRealFunc 	solFixedPolyInterp0 = sol.getSolAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solFixedPolyInterp1 = sol.getSolAsPolyInterp(1, 3);
	
	Visualizer::VisualizeRealFunction(solFixedPolyInterp0, "Forced Pendulum - angle in time",
																		t1, t2, 1000, "forced_pendulum_angle.txt");
	// shown together
	Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{ &solFixedPolyInterp0, &solFixedPolyInterp1 },
																				"Forced Pendulum - both variables", { "Angle", "Ang.vel." },
																				t1, t2, 1000, 
																				"forced_pendulum_multi_real_func.txt");

	Matrix<Real> curve_points(t_vals.size(), 2);
	for (int i = 0; i < t_vals.size(); i++)
	{
		curve_points(i, 0) = y1_fixed[i];
		curve_points(i, 1) = y2_fixed[i];
	}

	SplineInterpParametricCurve<2> phaseSpaceTrajectory(0.0, 1.0, curve_points);

	Visualizer::VisualizeParamCurve2D(phaseSpaceTrajectory, "Forced pendulum - phase space trajectory",
																		0.0, 1.0, t_vals.size(), "forced_pendulum_phase_space.txt");
}

void Demo_Damped_forced_pendulum()
{
	DampedForcedPendulumODE dampedForcedPen(9.8, 0.5, 1.2, 2.0/3);

	Real	t1 = 0.0, t2 = 10.0;
	Real	initAngle = 0.2;
	Vector<Real>	initCond{ initAngle, 0.0 };

	ODESystemSolver<RK5_CashKarp_Stepper> odeSolver(dampedForcedPen);
	ODESystemSolution sol = odeSolver.integrate(initCond, t1, t2, 0.01, 1e-09, 0.001);

	Vector<Real> t_vals{ sol.getTValues() }, y1{ sol.getXValues(0) }, y2{ sol.getXValues(1) };

	PolynomInterpRealFunc 	solY1 = sol.getSolAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solY2 = sol.getSolAsPolyInterp(1, 3);

	Visualizer::VisualizeRealFunction(solY1, "Damped Forced Pendulum - angle in time",
																		t1, t2, 1000, "damped_forced_pendulum_angle.txt");

	// shown together
	Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{ &solY1, &solY2 },
																				"Damped Forced Pendulum - both variables", { "Angle", "Ang.vel." },
																				t1, t2, 1000,
																				"damped_forced_pendulum_multi_real_func.txt");

	// visualize phase space trajectory as ParametricCurve2D
	Matrix<Real> curve_points(t_vals.size(), 2);
	for (int i = 0; i < t_vals.size(); i++)
	{
		curve_points(i, 0) = y1[i];
		curve_points(i, 1) = y2[i];
	}
	SplineInterpParametricCurve<2> phaseSpaceTrajectory(0.0, 1.0, curve_points);
	Visualizer::VisualizeParamCurve2D(phaseSpaceTrajectory, "Damped forced pendulum - phase space trajectory",
																		0.0, 1.0, 1000, "damped_forced_pendulum_phase_space.txt");
}

void Demo_Damped_forced_pendulum_compare_close_init_cond()
{
	DampedForcedPendulumODE dampedForcedPen(1.0, 0.1, 0.5, 2.0 / 3);

	Real	t1 = 0.0, t2 = 10.0;
	Real	initAngle1 = 0.2;
	Real	initAngle2 = 0.21;
	Vector<Real>	initCond1{ initAngle1, 0.0 };
	Vector<Real>	initCond2{ initAngle2, 0.0 };

	ODESystemSolver<RK5_CashKarp_Stepper> odeSolver(dampedForcedPen);
	
	ODESystemSolution sol1 = odeSolver.integrate(initCond1, t1, t2, 0.01, 1e-09, 0.001);
	ODESystemSolution sol2 = odeSolver.integrate(initCond2, t1, t2, 0.01, 1e-09, 0.001);

	std::cout << "Damped forced pendulum - comparing two close initial conditions" << std::endl;

	Vector<Real> t_vals1{ sol1.getTValues() }, y11{ sol1.getXValues(0) }, y21{ sol1.getXValues(1) };
	Vector<Real> t_vals2{ sol2.getTValues() }, y12{ sol2.getXValues(0) }, y22{ sol2.getXValues(1) };

	PolynomInterpRealFunc 	solY1 = sol1.getSolAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solW1 = sol1.getSolAsPolyInterp(1, 3);
	PolynomInterpRealFunc 	solY2 = sol2.getSolAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solW2 = sol2.getSolAsPolyInterp(1, 3);

	// shown together
	Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{ &solY1, & solY2 },
		"Damped Forced Pendulum - both cases", { "Angle1", "Angle2" },
		t1, t2, 1000,
		"damped_forced_pendulum_multi_real_func.txt");

	// visualize both phase space trajectories as ParametricCurve2D
	Matrix<Real> curve_points1(t_vals1.size(), 2);
	Matrix<Real> curve_points2(t_vals2.size(), 2);
	for (int i = 0; i < t_vals1.size(); i++)
	{
		curve_points1(i, 0) = y11[i];
		curve_points1(i, 1) = y21[i];
	}
	for (int i = 0; i < t_vals2.size(); i++)
	{
		curve_points2(i, 0) = y12[i];
		curve_points2(i, 1) = y22[i];
	}
	SplineInterpParametricCurve<2> phaseSpaceTrajectory1(0.0, 1.0, curve_points1);
	SplineInterpParametricCurve<2> phaseSpaceTrajectory2(0.0, 1.0, curve_points2);

	// visualize both phase space trajectories
	Visualizer::VisualizeMultiParamCurve2D(std::vector<IRealToVectorFunction<2>*>{ &phaseSpaceTrajectory1, &phaseSpaceTrajectory2 },
																				"Damped forced pendulum - phase space trajectories",
																				0.0, 1.0, 1000,
																				"damped_forced_pendulum_phase_space_compare.txt");
}
