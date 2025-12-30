#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "mml/MMLBase.h"

#include "mml/core/Derivation.h"

#include "mml/algorithms/ODESystemSolver.h"
#include "mml/algorithms/ODESystemStepCalculators.h"
#include "mml/algorithms/ODESystemSteppers.h"

#include "mml/tools/ConsolePrinter.h"
#include "mml/tools/Visualizer.h"
#include "mml/tools/Serializer.h"

#include "mpl/Mechanics/Pendulums.h"
#endif

#include <iostream>

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
	std::vector<ColumnFormat>				vecNames{ ColumnFormat("t", 11, 2, 'F'),
																				ColumnFormat("theta 1", 15, 8, 'F'),
																				ColumnFormat("theta 2", 15, 8, 'F'), };
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
