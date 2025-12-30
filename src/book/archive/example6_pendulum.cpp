#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/BaseUtils.h"
#include "base/Vector.h"

#include "tools/ConsolePrinter.h"
#include "tools/Visualizer.h"
#include "tools/Serializer.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/ODESystemSteppers.h"

#include "algorithms/RootFinding.h"
#include "algorithms/Statistics.h"
#include "algorithms/FunctionsAnalyzer.h"

#include "mpl/Mechanics/Pendulums.h"
#include "mpl/Mechanics/HarmonicOscillator.h"
#include "mpl/Mechanics/Oscillators.h"
#endif

using namespace MML;
using namespace MML::RootFinding;

using namespace MPL;

void Demo_SimplePendulum_Euler()
{
	// alternative way of defining the system with lambda function
	ODESystem pendSys(2, [](Real t, const Vector<Real>& x, Vector<Real>& dxdt)
		{
			Real Length = 1.0;
			dxdt[0] = x[1];
			dxdt[1] = -9.81 / Length * sin(x[0]);
		}
	);

	// setting parameters for our simulation
	Real	t1 = 0.0, t2 = 10.0;
	int   expectNumSteps = 1000;							// means our step = 10 / 1000 = 0.001
	Real	initAngle = Utils::DegToRad(30);		// initial angle set to 30 degrees
	Vector<Real>	initCond{ initAngle, 0.0 };

	// solving and getting solutions
	ODESystemFixedStepSolver	fixedSolver(pendSys, StepCalculators::EulerStepCalc);
	ODESystemSolution sol = fixedSolver.integrate(initCond, t1, t2, expectNumSteps);

	Vector<Real> sol_x{ sol.getTValues() }, sol_y1{ sol.getXValues(0) }, sol_y2{ sol.getXValues(1) };

	// console visualization
	std::cout << "\n\n****  Euler method - fixed stepsize  **********\n";
	std::vector<ColDesc> vecNames{ ColDesc("t", 11, 2, 'F'), ColDesc("angle", 15, 8, 'F'), ColDesc("ang.vel.", 15, 8, 'F') };
	std::vector<Vector<Real>*>	vecVals{ &sol_x, &sol_y1, &sol_y2 };

	VectorTablePrinter	vvp(vecNames, vecVals);

	vvpprint();

	// visualizing solutions
	PolynomInterpRealFunc 	solAngle    = sol.getSolAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solAngleVel = sol.getSolAsPolyInterp(1, 3);

	Visualizer::VisualizeRealFunction(solAngle, "Pendulum - angle in time",
																		t1, t2, 200, "pendulum_angle_euler.txt");

	Visualizer::VisualizeRealFunction(solAngleVel, "Pendulum - ang.vel. in time",
																		t1, t2, 200, "pendulum_ang_vel_euler.txt");

	// shown together as multi-real function
	Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{ &solAngle, &solAngleVel },
																				"Pendulum - both variables", { "Angle", "Ang.vel." },
																				t1, t2, 200,
																				"pendulum_multi_real_func_euler.txt");

	// visualized together as ODE system solution
	Visualizer::VisualizeODESysSolAsMultiFunc(sol, "Pendulum - Euler method", std::vector<std::string>{"angle", "angle.vel."},
																						"pendulum_euler.txt");
}

void ComparePendulumPeriods(
	PendulumODE& sys,
	Vector<Real>& initCond,
	Real t1, Real t2, int expectNumSteps, Real minSaveInterval)
{
	ODESystemFixedStepSolver	fixedSolver(sys, StepCalculators::RK4_Basic);
	ODESystemSolver<RK5_CashKarp_Stepper> adaptSolver(sys);

	Vector<Real> initAnglesDeg{ 1, 2, 5, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0 };
	Vector<Real> initAngles(initAnglesDeg.size());
	for (int i = 0; i < initAngles.size(); i++)
		initAngles[i] = Utils::DegToRad(initAnglesDeg[i]);

	TablePrinter<double, double> data("Angle", 5, 1,
		{ "Linear ", "Exact  ", "Diff %", "Fix.step.sim", "Diff %", "Adapt.step.sim", "Diff % " },
		{ {10,6,'F'}, {10,6,'F'}, {7,3,'F'}, {13,7,'F'}, {8,5,'F'}, {13,7,'F'}, {8,5,'F'} });

	for (int i = 0; i < initAngles.size(); i++)
	{
		initCond[0] = initAngles[i];

		// solving with fixed and adaptive solver
		ODESystemSolution solF = fixedSolver.integrate(initCond, t1, t2, expectNumSteps);
		ODESystemSolution solA = adaptSolver.integrate(initCond, t1, t2, minSaveInterval, 1e-06, 0.01);

		// calculating period from solutions
		PolynomInterpRealFunc solFInterp = solF.getSolAsPolyInterp(0, 3);
		PolynomInterpRealFunc solAInterp = solA.getSolAsPolyInterp(0, 3);
		RealFunctionAnalyzer fa3(solFInterp);
		RealFunctionAnalyzer fa4(solAInterp);
		Real periodSimulFixed = fa3.calcRootsPeriod(t1, t2, 1000);
		Real periodSimulAdapt = fa4.calcRootsPeriod(t1, t2, 1000);

		// analytical formulas
		Real periodLin = sys.calcPeriodLinearized();
		Real periodExact = sys.calcExactPeriod(initAngles[i]);
		Real diffPercent = Abs(periodLin - periodExact) / periodExact * 100;
		Real diffPercentFixed = Abs(2 * periodSimulFixed - periodExact) / periodExact * 100;
		Real diffPercentAdapt = Abs(2 * periodSimulAdapt - periodExact) / periodExact * 100;

		data.addRow(initAnglesDeg[i], { double(periodLin), double(periodExact), double(diffPercent),
				double(2 * periodSimulFixed), double(diffPercentFixed), double(2 * periodSimulAdapt), double(diffPercentAdapt) });
	}

	dataprint();
}

void CompareAdaptiveStepCounts(
	PendulumODE& sys,
	Vector<Real> initCond,
	Real t1, Real t2, Real minSaveInterval)
{
	Vector<Real> acc{ 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10 };
	Vector<Real> numSteps(acc.size());

	std::cout << "\nAccuracy   Num steps   OK steps   Bad steps" << std::endl;
	for (int i = 0; i < acc.size(); i++)
	{
		ODESystemSolver<RK5_CashKarp_Stepper> adaptSolver(sys);
		ODESystemSolution solAdapt = adaptSolver.integrate(initCond, t1, t2, minSaveInterval, acc[i], 0.01);
		numSteps[i] = solAdapt.getTotalNumSteps();

		std::cout << std::setw(7) << acc[i] << "        " << std::setw(3) << numSteps[i] << "        "
			<< std::setw(3) << solAdapt.getNumStepsOK() << "        "
			<< std::setw(3) << solAdapt.getNumStepsBad() << std::endl;
	}
}

void CompareFixedAndAdaptiveSolver(PendulumODE& sys,
																		Vector<Real> initCond,
																		Real t1, Real t2, int expectNumSteps)
{
	Real	minSaveInterval = (t2 - t1) / expectNumSteps;

	ODESystemFixedStepSolver	fixedSolver(sys, StepCalculators::RK4_Basic);
	ODESystemSolver<RK5_CashKarp_Stepper> adaptSolver(sys);

	ODESystemSolution solFixed = fixedSolver.integrate(initCond, t1, t2, expectNumSteps);
	ODESystemSolution solAdapt = adaptSolver.integrate(initCond, t1, t2, minSaveInterval, 1e-08, 0.01);

	Vector<Real> x_fixed{ solFixed.getTValues() }, y1_fixed{ solFixed.getXValues(0) }, y2_fixed{ solFixed.getXValues(1) };
	Vector<Real> x_adapt{ solAdapt.getTValues() }, y1_adapt{ solAdapt.getXValues(0) }, y2_adapt{ solAdapt.getXValues(1) };
	Vector<Real> y1_adapt_bounded_pi_pi, y1_adapt_bounded_0_2pi;	// solution with angle reduced to range {-Pi, Pi} and {0, 2pi}

	for (int i = 0; i < solAdapt.getXValues(0).size(); i++)
	{
		Real angle = solAdapt.getXValues(0)[i];

		y1_adapt_bounded_pi_pi.push_back(Utils::AngleToPiPiRange(angle));
		y1_adapt_bounded_0_2pi.push_back(Utils::AngleTo2PiRange(angle));
	}

	std::cout << "\n\n****  Runge-Kutta 4th order - fixed stepsize  **********  Runge-Kutta 4th order - adaptive stepper  ****\n";
	std::vector<ColDesc>				vecNames{ ColDesc("t", 11, 2, 'F'), ColDesc("angle", 15, 8, 'F'), ColDesc("ang.vel.", 15, 8, 'F'),
																				ColDesc("t", 22, 2, 'F'), ColDesc("angle", 12, 8, 'F'), ColDesc("ang.vel.", 12, 8, 'F') };
	std::vector<Vector<Real>*>	vecVals{ &x_fixed, &y1_fixed, &y2_fixed,
																			 &x_adapt, &y1_adapt, &y2_adapt };
	VectorTablePrinter	vvp(vecNames, vecVals);

	//vvpprint();

	// getting solutions as polynomials
	PolynomInterpRealFunc 	angleFixedSol = solFixed.getSolAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	angleVelFixedSol = solFixed.getSolAsPolyInterp(1, 3);

	PolynomInterpRealFunc 	angleAdaptSol = solAdapt.getSolAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	angleVelAdaptSol = solAdapt.getSolAsPolyInterp(1, 3);
	LinearInterpRealFunc		angleAdaptBounded(solAdapt.getTValues(), y1_adapt_bounded_pi_pi);

	Visualizer::VisualizeRealFunction(solAdapt.getSolAsPolyInterp(0, 3), "Pendulum - angle in time",
		t1, t2, 500, "pendulum_angle.txt");

	Visualizer::VisualizeRealFunction(angleAdaptBounded, "Pendulum - bounded angle in time",
		t1, t2, 500, "pendulum_bounded_angle.txt");

	// shown together
	Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{ &angleAdaptSol, & angleVelAdaptSol },
																				"Pendulum - both variables", { "Angle", "Ang.vel" },
																				t1, t2, 200,
																				"pendulum_multi_real_func.txt");


	LinearInterpRealFunc angleVelFunc(solAdapt.getTValues(), y2_adapt);

	// show together
	Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{ &angleAdaptBounded, & angleVelFunc },
																				"Pendulum - both variables (bounded angle)", { "Angle", "Ang.vel." },
																				t1, t2, 500,
																				"pendulum_multi_real_func_bounded_angle.txt");

	Matrix<Real> curve_points(x_adapt.size(), 2);
	for (int i = 0; i < x_adapt.size(); i++)
	{
		curve_points(i, 0) = y1_adapt[i];
		curve_points(i, 1) = y2_adapt[i];
	}
	SplineInterpParametricCurve<2> phaseSpaceTrajectory(0.0, 1.0, curve_points);
	Visualizer::VisualizeParamCurve2D(phaseSpaceTrajectory, "Pendulum - phase space trajectory",
																		0.0, 1.0, 1000, "pendulum_phase_space.txt");

	Real periodLinear = sys.calcPeriodLinearized();
	Real exactPeriod = sys.calcExactPeriod(initCond[0]);
	// extracting period T from solutions
	Real simulPeriodFixed = RealFunctionAnalyzer(angleFixedSol).calcRootsPeriod(t1, t2, 1000);
	Real simulPeriodAdapt = RealFunctionAnalyzer(angleAdaptSol).calcRootsPeriod(t1, t2, 1000);

	std::cout << "Pendulum period approx. linear : " << periodLinear << std::endl;
	std::cout << "Pendulum period analytic exact : " << exactPeriod << std::endl;
	std::cout << "Pendulum period RK4 fixed step : " << 2 * simulPeriodFixed << std::endl;
	std::cout << "Pendulum period RK4 adapt.step : " << 2 * simulPeriodAdapt << std::endl;
}

void ComparePhaseSpace(PendulumODE& sys,
												Vector<Real> initCond,
												Real t1, Real t2, Real minSaveInterval)
{
	ODESystemFixedStepSolver	fixedSolver(sys, StepCalculators::RK4_Basic);
	ODESystemSolver<RK5_CashKarp_Stepper> adaptSolver(sys);

	Real initVel = -2.0;

	Vector<Real> initAnglesDeg{ 10.0, 30.0, 50.0, 70.0, 90.0, 120.0, 150, 175 };
	Vector<Real> initAngles(initAnglesDeg.size());
	for (int i = 0; i < initAngles.size(); i++)
		initAngles[i] = Utils::DegToRad(initAnglesDeg[i]);

	std::vector<IRealToVectorFunction<2> *> vecPhaseSpaceTraj;
	std::vector<std::string> legends;

	for (int i = 0; i < initAngles.size(); i++)
	{
		initCond[0] = initAngles[i];
		initCond[1] = initVel;

		// check if it is a closed orbit
		bool isClosedOrbit = sys.isClosedOrbit(initCond[0], initCond[1]);
		std::cout << "Pendulum with initial angle " << initAnglesDeg[i] << " degrees and initial velocity "
			<< initVel << " is closed orbit: " << (isClosedOrbit ? "yes" : "no") << std::endl;

		ODESystemSolution solAdapt = adaptSolver.integrate(initCond, t1, t2, minSaveInterval, 1e-08, 0.01);

		Vector<Real> x_adapt{ solAdapt.getTValues() }, y1_adapt{ solAdapt.getXValues(0) }, y2_adapt{ solAdapt.getXValues(1) };
		Vector<Real> y1_adapt_bounded_pi_pi, y1_adapt_bounded_0_2pi;	// solution with angle reduced to range {-Pi, Pi} and {0, 2pi}

		for (int j = 0; j < solAdapt.getXValues(0).size(); j++)
		{
			Real angle = solAdapt.getXValues(0)[j];

			y1_adapt_bounded_pi_pi.push_back(Utils::AngleToPiPiRange(angle));
			y1_adapt_bounded_0_2pi.push_back(Utils::AngleTo2PiRange(angle));
		}

		// creating phase space trajectory as ParametricCurve2D
		Matrix<Real> curve_points_bounded_pi_pi;
		if (isClosedOrbit)
		{
			Matrix<Real> curve_points(x_adapt.size(), 2);
			curve_points_bounded_pi_pi.Resize(x_adapt.size(), 2);

			for (int i = 0; i < x_adapt.size(); i++)
			{
				curve_points(i, 0) = y1_adapt[i];
				curve_points(i, 1) = y2_adapt[i];

				curve_points_bounded_pi_pi(i, 0) = y1_adapt_bounded_pi_pi[i];
				//			curve_points2(i, 0) = y1_adapt_bounded_0_2pi[i];
				curve_points_bounded_pi_pi(i, 1) = y2_adapt[i];
			}
		}
		else
		{
			// BUT, if it is NOT a closed orbit, we need to extract one single period
			// first, and find first point where we go from -pi to pi
			int firstCrossingIndex = -1;
			for (int j = 0; j < y1_adapt_bounded_pi_pi.size() - 1; j++)
			{
				if (y1_adapt_bounded_pi_pi[j] < -Constants::PI/2 && y1_adapt_bounded_pi_pi[j + 1] > Constants::PI/2)
				{
					firstCrossingIndex = j+1;
					break;
				}
			}
			// now second point
			int secondCrossingIndex = -1;
			for (int j = firstCrossingIndex+1; j < y1_adapt_bounded_pi_pi.size() - 1; j++)
			{
				if (y1_adapt_bounded_pi_pi[j] < -Constants::PI/2 && y1_adapt_bounded_pi_pi[j + 1] > Constants::PI/2)
				{
					secondCrossingIndex = j;
					break;
				}
			}
			// now we need to extract only those points to our curve_points_bounded_pi_pi
			if (firstCrossingIndex < 0 || secondCrossingIndex < 0)
			{
				std::cerr << "Error: could not find crossing points for non-closed orbit!" << std::endl;
				continue;
			}
			int numPoints = secondCrossingIndex - firstCrossingIndex + 1;
			curve_points_bounded_pi_pi.Resize(numPoints, 2);
			for (int j = 0; j < numPoints; j++)
			{
				curve_points_bounded_pi_pi(j, 0) = y1_adapt_bounded_pi_pi[firstCrossingIndex + j];
				curve_points_bounded_pi_pi(j, 1) = y2_adapt[firstCrossingIndex + j];
			}
		}

		LinInterpParametricCurve<2>* phaseSpaceTrajectory = new LinInterpParametricCurve<2>(curve_points_bounded_pi_pi);

		vecPhaseSpaceTraj.push_back(phaseSpaceTrajectory);
		legends.push_back("phi0 = " + std::to_string((int) initAnglesDeg[i]) + "; w0 = " + std::to_string(initVel));
	}

	// visualizing phase space trajectories as ParametricCurve2D
	Visualizer::VisualizeMultiParamCurve2D(vecPhaseSpaceTraj, legends,
																				0.0, 1.0, 1000, "pendulum_phase_space_multi.txt");
}

void Demo_SimplePendulum()
{
	ODESystem pendSys(2, [](Real t, const Vector<Real>& x, Vector<Real>& dxdt)
		{
			Real Length = 1.0;
			dxdt[0] = x[1];
			dxdt[1] = -9.81 / Length * sin(x[0]);
		}
	);

	Real	pendulumLen = 1.0;
	PendulumODE		sys = PendulumODE(pendulumLen);

	Real	t1 = 0.0, t2 = 5.0;
	int   expectNumSteps = 1000;
	Real	minSaveInterval = (t2 - t1) / expectNumSteps;
	
	Real	initAngle = Utils::DegToRad(170), initAngVel = -4.5;
	Vector<Real>	initCond{ initAngle, initAngVel };

	//CompareFixedAndAdaptiveSolver(sys, initCond, t1, t2, expectNumSteps);

	// Comparing periods for different initial angles
	ComparePendulumPeriods(sys, initCond, t1, t2, expectNumSteps, minSaveInterval);

	// Comparing number of adaptive steps needed, for given accuracy
	CompareAdaptiveStepCounts(sys, initCond, t1, t2, minSaveInterval);

	ComparePhaseSpace(sys, initCond, t1, t2, minSaveInterval);
}

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

void Demo_HarmonicOscilator()
{
	HarmonicOscillator hoSys(1.0, 0.5);

	Real	t1 = 0.0, t2 = 30.0;
	Real  initVel = 0.0;
	Real	initAngle = 0.5;

	// using functions provided by HarmonicOscillator class
	RealFunctionFromStdFunc y1 = hoSys.getPositionFunction(initAngle, initVel);
	RealFunctionFromStdFunc y2 = hoSys.getVelocityFunction(initAngle, initVel);

	Visualizer::VisualizeRealFunction(y1, "Harmonic Oscilator - angle in time",
																		t1, t2, 200, "harmonic_oscillator_angle.txt");
	// shown together
	Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{ &y1, &y2 },
																				"Harmonic Oscilator - both variables", { "Angle", "Ang.vel." },
																				t1, t2, 200,
																				"harmonic_oscillator_multi_real_func.txt");
	
	// now we have to visualize phase space trajectory as ParametricCurve2D
	int numCurvePoints = 1001;	// number of points in phase space trajectory
	Matrix<Real> curve_points(numCurvePoints, 2);	// 1000 points for phase space trajectory
	for (int i = 0; i < numCurvePoints; i++)
	{
		Real t = t1 + (t2 - t1) * i / (numCurvePoints - 1);
		curve_points(i, 0) = y1(t);	// position
		curve_points(i, 1) = y2(t);	// velocity
	}
	SplineInterpParametricCurve<2> phaseSpaceTrajectory(0.0, 1.0, curve_points);

	Visualizer::VisualizeParamCurve2D(phaseSpaceTrajectory, "Harmonic Oscilator, exact - phase space trajectory",
																		0.0, 1.0, numCurvePoints, "harmonic_oscillator_phase_space_exact.txt");
	

	// solving differential equation
	int   expectNumSteps = 1000;
	Vector<Real>	initCond{ initAngle, initVel };

	ODESystemFixedStepSolver	fixedSolver(hoSys, StepCalculators::RK4_Basic);
	ODESystemSolution sol = fixedSolver.integrate(initCond, t1, t2, expectNumSteps);

	Vector<Real> t_vals{ sol.getTValues() }, y1_fixed{ sol.getXValues(0) }, y2_fixed{ sol.getXValues(1) };

	PolynomInterpRealFunc 	solInterpY1 = sol.getSolAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solInterpY2 = sol.getSolAsPolyInterp(1, 3);

	Matrix<Real> curve_points2(t_vals.size(), 2);
	for (int i = 0; i < t_vals.size(); i++)
	{
		curve_points2(i, 0) = y1_fixed[i];
		curve_points2(i, 1) = y2_fixed[i];
	}
	SplineInterpParametricCurve<2> phaseSpaceTrajectory2(0.0, 1.0, curve_points2);
	Visualizer::VisualizeParamCurve2D(phaseSpaceTrajectory2, "Harmonic Oscilator, diff.eq.sol. - phase space trajectory",
																		0.0, 1.0, t_vals.size(), "harmonic_oscillator_phase_space_diff_eq_soluton.txt");
}

void Demo_DampedHarmonicOscilator()
{
	DampedHarmonicOscillator hoSys(1.0, 0.5, 0.8);
	Real	t1 = 0.0, t2 = 30.0;
	Real  initVel = 0.0;
	Real	initAngle = 0.5;

	// using functions provided by HarmonicOscillator class
	RealFunctionFromStdFunc y1 = hoSys.getPositionFunction(initAngle, initVel);
	RealFunctionFromStdFunc y2 = hoSys.getVelocityFunction(initAngle, initVel);
	
	Visualizer::VisualizeRealFunction(y1, "Damped Harmonic Oscilator - angle in time",
																		t1, t2, 200, "damped_harmonic_oscillator_angle.txt");
	// shown together
	Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{ &y1, &y2 },
																				"Damped Harmonic Oscilator - both variables", { "Angle", "Ang.vel." },
																				t1, t2, 200,
																				"damped_harmonic_oscillator_multi_real_func.txt");

	// now we have to visualize phase space trajectory as ParametricCurve2D
	int numCurvePoints = 1001;	// number of points in phase space trajectory
	Matrix<Real> curve_points(numCurvePoints, 2);	// 1000 points for phase space trajectory
	for (int i = 0; i < numCurvePoints; i++)
	{
		Real t = t1 + (t2 - t1) * i / (numCurvePoints - 1);
		curve_points(i, 0) = y1(t);	// position
		curve_points(i, 1) = y2(t);	// velocity
	}

	SplineInterpParametricCurve<2> phaseSpaceTrajectory(0.0, 1.0, curve_points);
	Visualizer::VisualizeParamCurve2D(phaseSpaceTrajectory, "Damped Harmonic Oscilator - phase space trajectory",
																		0.0, 1.0, numCurvePoints, "damped_harmonic_oscillator_phase_space_exact.txt");


	// solving differential equation
	int   expectNumSteps = 1000;
	Vector<Real>	initCond{ initAngle, initVel };

	ODESystemFixedStepSolver	fixedSolver(hoSys, StepCalculators::RK4_Basic);
	ODESystemSolution sol = fixedSolver.integrate(initCond, t1, t2, expectNumSteps);

	Vector<Real> t_vals{ sol.getTValues() }, y1_fixed{ sol.getXValues(0) }, y2_fixed{ sol.getXValues(1) };

	PolynomInterpRealFunc 	solInterpY1 = sol.getSolAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solInterpY2 = sol.getSolAsPolyInterp(1, 3);
	
	//Visualizer::VisualizeRealFunction(solInterpY1, "Damped Harmonic Oscilator - angle in time",
	//																	t1, t2, 500, "damped_harmonic_oscillator_angle_diff_eq.txt");

	//// shown together
	//Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{ &solInterpY1, &solInterpY2 },
	//																			"Damped Harmonic Oscilator - both variables", { "Angle", "Ang.vel." },
	//																			t1, t2, 500,
	//																			"damped_harmonic_oscillator_multi_real_func_diff_eq.txt");

	// visualize phase space trajectory as ParametricCurve2D
	Matrix<Real> curve_points2(t_vals.size(), 2);
	for (int i = 0; i < t_vals.size(); i++)
	{
		curve_points2(i, 0) = y1_fixed[i];
		curve_points2(i, 1) = y2_fixed[i];
	}
	SplineInterpParametricCurve<2> phaseSpaceTrajectory2(0.0, 1.0, curve_points2);
	Visualizer::VisualizeParamCurve2D(phaseSpaceTrajectory2, "Damped Harmonic Oscilator - phase space trajectory",
																		0.0, 1.0, t_vals.size(), "damped_harmonic_oscillator_phase_space_diff_eq.txt");
}

void Demo_ForcedHarmonicOscilator()
{
	ForcedHarmonicOscillator hoSys(1.0, 0.5, 0.1, 0.707);

	std::cout << "Forced Harmonic Oscilator - using functions provided by class" << std::endl;

	std::cout << "Natural frequency : " << hoSys.getNaturalFrequency() << std::endl;
	std::cout << "Resonance frequency : " << hoSys.getResonanceFrequency() << std::endl;
	std::cout << "Period : " << hoSys.getPeriod() << std::endl;

	Real	t1 = 0.0, t2 = 100.0;
	Real  initVel = 0.0;
	Real	initAngle = 0.5;

	// using functions provided by HarmonicOscillator class
	RealFunctionFromStdFunc y1 = hoSys.getPositionFunction(initAngle, initVel);
	RealFunctionFromStdFunc y2 = hoSys.getVelocityFunction(initAngle, initVel);
	Visualizer::VisualizeRealFunction(y1, "Forced Harmonic Oscilator - angle in time",
																		t1, t2, 200, "forced_harmonic_oscillator_angle.txt");
	// shown together
	Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{ &y1, &y2 },
																				"Forced Harmonic Oscilator - both variables", { "Angle", "Ang.vel." },
																				t1, t2, 200,
																				"forced_harmonic_oscillator_multi_real_func.txt");
	
	// now we have to visualize phase space trajectory as ParametricCurve2D
	int numCurvePoints = 1001;	// number of points in phase space trajectory
	Matrix<Real> curve_points(numCurvePoints, 2);	// 1000 points for phase space trajectory
	for (int i = 0; i < numCurvePoints; i++)
	{
		Real t = t1 + (t2 - t1) * i / (numCurvePoints - 1);
		curve_points(i, 0) = y1(t);	// position
		curve_points(i, 1) = y2(t);	// velocity
	}
	SplineInterpParametricCurve<2> phaseSpaceTrajectory(0.0, 1.0, curve_points);
	Visualizer::VisualizeParamCurve2D(phaseSpaceTrajectory, "Forced Harmonic Oscilator - phase space trajectory",
																		0.0, 1.0, numCurvePoints, "forced_harmonic_oscillator_phase_space_exact.txt");


	// now we will solve differential equation
	Real	t3 = 0.0, t4 = 100.0;
	int   expectNumSteps = 1000;
	Vector<Real>	initCond{ initAngle, initVel };
	ODESystemFixedStepSolver	fixedSolver(hoSys, StepCalculators::RK4_Basic);
	ODESystemSolution sol = fixedSolver.integrate(initCond, t3, t4, expectNumSteps);
	Vector<Real> t_vals{ sol.getTValues() }, y1_fixed{ sol.getXValues(0) }, y2_fixed{ sol.getXValues(1) };
	PolynomInterpRealFunc 	solInterpY1 = sol.getSolAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solInterpY2 = sol.getSolAsPolyInterp(1, 3);

	Visualizer::VisualizeRealFunction(solInterpY1, "Forced Harmonic Oscilator - angle in time",
		t3, t4, 200, "forced_harmonic_oscillator_angle_diff_eq.txt");

	// shown together
	Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{ &solInterpY1, &solInterpY2 },
																				"Forced Harmonic Oscilator - both variables", { "Angle", "Ang.vel." },
																				t3, t4, 200,
																				"forced_harmonic_oscillator_multi_real_func_diff_eq.txt");

	// visualize phase space trajectory as ParametricCurve2D
	Matrix<Real> curve_points2(t_vals.size(), 2);
	for (int i = 0; i < t_vals.size(); i++)
	{
		curve_points2(i, 0) = y1_fixed[i];
		curve_points2(i, 1) = y2_fixed[i];
	}
	SplineInterpParametricCurve<2> phaseSpaceTrajectory2(0.0, 1.0, curve_points2);
	Visualizer::VisualizeParamCurve2D(phaseSpaceTrajectory2, "Forced Harmonic Oscilator - phase space trajectory",
																		0.0, 1.0, t_vals.size(), "forced_harmonic_oscillator_phase_space_diff_eq.txt");
}

void Demo_DampedForcedHarmonicOscilator()
{
	DampedForcedHarmonicOscillator hoSys(1.0, 0.5, 0.1, 0.1, 0.95);

	std::cout << "Damped Forced Harmonic Oscilator - using functions provided by class" << std::endl;
	std::cout << "Natural frequency : " << hoSys.getNaturalFrequency() << std::endl;
	std::cout << "Period : " << hoSys.getPeriod() << std::endl;
	std::cout << "Damping ratio : " << hoSys.getDampingRatio() << std::endl;	

	Real	t1 = 0.0, t2 = 30.0;
	Real  initVel = 0.0;
	Real	initAngle = 0.5;
	
	// using functions provided by HarmonicOscillator class
	RealFunctionFromStdFunc y1 = hoSys.getPositionFunction(initAngle, initVel);
	RealFunctionFromStdFunc y2 = hoSys.getVelocityFunction(initAngle, initVel);
	
	Visualizer::VisualizeRealFunction(y1, "Damped Forced Harmonic Oscilator - angle in time",
																		t1, t2, 200, "damped_forced_harmonic_oscillator_angle.txt");
	
	// shown together
	Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{ &y1, &y2 },
																				"Damped Forced Harmonic Oscilator - both variables", { "Angle", "Ang.vel." },
																				t1, t2, 200,
																				"damped_forced_harmonic_oscillator_multi_real_func.txt");
	
	// now we have to visualize phase space trajectory as ParametricCurve2D
	int numCurvePoints = 1001;	// number of points in phase space trajectory
	Matrix<Real> curve_points(numCurvePoints, 2);	// 1000 points for phase space trajectory
	for (int i = 0; i < numCurvePoints; i++)
	{
		Real t = t1 + (t2 - t1) * i / (numCurvePoints - 1);
		curve_points(i, 0) = y1(t);	// position
		curve_points(i, 1) = y2(t);	// velocity
	}

	SplineInterpParametricCurve<2> phaseSpaceTrajectory(0.0, 1.0, curve_points);
	
	Visualizer::VisualizeParamCurve2D(phaseSpaceTrajectory, "Damped Forced Harmonic Oscilator - phase space trajectory",
																		0.0, 1.0, numCurvePoints, "damped_forced_harmonic_oscillator_phase_space_exact.txt");

	// now we will solve differential equation
	int   expectNumSteps = 2000;
	Vector<Real>	initCond{ initAngle, initVel };
	
	ODESystemFixedStepSolver	fixedSolver(hoSys, StepCalculators::RK4_Basic);
	ODESystemSolution sol = fixedSolver.integrate(initCond, t1, t2, expectNumSteps);
	
	Vector<Real> t_vals{ sol.getTValues() }, y1_fixed{ sol.getXValues(0) }, y2_fixed{ sol.getXValues(1) };
	
	PolynomInterpRealFunc 	solInterpY1 = sol.getSolAsPolyInterp(0, 3);
	PolynomInterpRealFunc 	solInterpY2 = sol.getSolAsPolyInterp(1, 3);

	Visualizer::VisualizeRealFunction(solInterpY1, "Damped Forced Harmonic Oscilator - angle in time",
																		t1, t2, 500, "damped_forced_harmonic_oscillator_angle_diff_eq.txt");

	// shown together
	Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{ &solInterpY1, &solInterpY2 },
																				"Damped Forced Harmonic Oscilator - both variables", { "Angle", "Ang.vel." },
																				t1, t2, 500,
																				"damped_forced_harmonic_oscillator_multi_real_func_diff_eq.txt");

	// visualize phase space trajectory as ParametricCurve2D
	Matrix<Real> curve_points2(t_vals.size(), 2);
	for (int i = 0; i < t_vals.size(); i++)
	{
		curve_points2(i, 0) = y1_fixed[i];
		curve_points2(i, 1) = y2_fixed[i];
	}
	SplineInterpParametricCurve<2> phaseSpaceTrajectory2(0.0, 1.0, curve_points2);
	
	Visualizer::VisualizeParamCurve2D(phaseSpaceTrajectory2, "Damped Forced Harmonic Oscilator - phase space trajectory",
																		0.0, 1.0, 500, "damped_forced_harmonic_oscillator_phase_space_diff_eq.txt");
}

void Example6_Pendulum()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                   EXAMPLE 5 - Pendulum                        ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	//Demo_SimplePendulum();
	//Demo_SimplePendulumEuler();
	//Demo_DampedPendulum();
	Demo_DampedPendulum_compare_three_cases();
	//Demo_ForcedPendulum();
	//Demo_Damped_forced_pendulum();
	//Demo_Damped_forced_pendulum_compare_close_init_cond();

	//Demo_HarmonicOscilator();
	//Demo_DampedHarmonicOscilator();
	//Demo_ForcedHarmonicOscilator();
	//Demo_DampedForcedHarmonicOscilator();
}