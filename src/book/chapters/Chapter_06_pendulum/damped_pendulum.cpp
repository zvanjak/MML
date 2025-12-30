#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "mml/MMLBase.h"

#include "mml/base/BaseUtils.h"
#include "mml/base/Vector.h"

#include "mml/tools/Visualizer.h"

#include "mml/algorithms/ODESystemSolver.h"
#include "mml/algorithms/ODESystemStepCalculators.h"
#include "mml/algorithms/ODESystemSteppers.h"

#include "mpl/Mechanics/Pendulums.h"
#endif

using namespace MML;
using namespace MPL;

void Demo_DampedPendulum()
{
	DampedPendulumODE odeSys(1.0, 0.5);

	Real	t1 = 0.0, t2 = 20.0;
	int   expectNumSteps = 1000;
	Real	initAngle = 0.5;
	Vector<Real>	initCond{ initAngle, 0.0 };
	
	ODESystemFixedStepSolver	fixedSolver(odeSys, StepCalculators::RK4_Basic);
	ODESystemSolution sol = fixedSolver.integrate(initCond, t1, t2, expectNumSteps);
	
	Vector<Real> t_vals{ sol.getTValues() }, y1_fixed{ sol.getXValues(0) }, y2_fixed{ sol.getXValues(1) };
	
	PolynomInterpRealFunc 	solInterpY1 = sol.getSolAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solInterpY2 = sol.getSolAsPolyInterp(1, 3);
	
	Visualizer::VisualizeRealFunction(solInterpY1, "Dumped Pendulum - angle in time",
																		t1, t2, 200, "dumped_pendulum_angle.txt");
	// shown together
	Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{ &solInterpY1, &solInterpY2 },
																					"Dumped Pendulum - both variables", { "Angle", "Ang.vel." },
																					t1, t2, 200, 
																					"dumped_pendulum_multi_real_func.txt");

	// visualize phase space trajectory as ParametricCurve2D
	Matrix<Real> curve_points(t_vals.size(), 2);
	for (int i = 0; i < t_vals.size(); i++)
	{
		curve_points(i, 0) = y1_fixed[i];
		curve_points(i, 1) = y2_fixed[i];
	}

	SplineInterpParametricCurve<2> phaseSpaceTrajectory(0.0, 1.0, curve_points);

	Visualizer::VisualizeParamCurve2D(phaseSpaceTrajectory, "Damped pendulum - phase space trajectory",
																		0.0, 1.0, t_vals.size(), "damped_pendulum_phase_space.txt");
}

void Demo_DampedPendulum_compare_three_cases()
{
	DampedPendulumODE testODESys(1.0, 0.5);

	Real criticalDamping = testODESys.calcCriticalDampingCoeff();

	DampedPendulumODE odeSys1(1.0, criticalDamping * 0.3);
	DampedPendulumODE odeSys11(1.0, criticalDamping * 0.5);
	DampedPendulumODE odeSys12(1.0, criticalDamping * 0.7);

	DampedPendulumODE odeSys2(1.0, criticalDamping);
	DampedPendulumODE odeSys3(1.0, criticalDamping * 1.2);

	Real	t1 = 0.0, t2 = 5;
	int   expectNumSteps = 1000;
	Real	initAngle = 0.5;
	Vector<Real>	initCond{ initAngle, 0.0 };

	ODESystemFixedStepSolver	fixedSolver1(odeSys1, StepCalculators::RK4_Basic);
	ODESystemFixedStepSolver	fixedSolver11(odeSys11, StepCalculators::RK4_Basic);
	ODESystemFixedStepSolver	fixedSolver12(odeSys12, StepCalculators::RK4_Basic);

	ODESystemFixedStepSolver	fixedSolver2(odeSys2, StepCalculators::RK4_Basic);
	ODESystemFixedStepSolver	fixedSolver3(odeSys3, StepCalculators::RK4_Basic);

	ODESystemSolution sol1  = fixedSolver1.integrate(initCond, t1, t2, expectNumSteps);
	ODESystemSolution sol11 = fixedSolver11.integrate(initCond, t1, t2, expectNumSteps);
	ODESystemSolution sol12 = fixedSolver12.integrate(initCond, t1, t2, expectNumSteps);

	ODESystemSolution sol2 = fixedSolver2.integrate(initCond, t1, t2, expectNumSteps);
	ODESystemSolution sol3 = fixedSolver3.integrate(initCond, t1, t2, expectNumSteps);

	PolynomInterpRealFunc 	solInterpY1_1  = sol1.getSolAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solInterpY1_11 = sol11.getSolAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solInterpY1_12 = sol12.getSolAsPolyInterp(0, 3);

	PolynomInterpRealFunc 	solInterpY1_2 = sol2.getSolAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solInterpY1_3 = sol3.getSolAsPolyInterp(0, 3);

	// show together
	Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{ &solInterpY1_1, & solInterpY1_11, & solInterpY1_12, &solInterpY1_2, &solInterpY1_3 },
																				"Damped Pendulum", { "Crit * 0.3", "Crit * 0.5", "Crit * 0.7", "Critical", "Over critical"},
																				t1, t2, 1000,
																				"damped_pendulum_compare_cases.txt");
}
