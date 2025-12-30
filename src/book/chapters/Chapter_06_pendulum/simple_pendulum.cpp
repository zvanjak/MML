#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "mml/MMLBase.h"

#include "mml/interfaces/IODESystem.h"

#include "mml/base/BaseUtils.h"
#include "mml/base/Vector.h"

#include "mml/tools/ConsolePrinter.h"
#include "mml/tools/Visualizer.h"

#include "mml/algorithms/ODESystemSolver.h"
#include "mml/algorithms/ODESystemStepCalculators.h"
#include "mml/algorithms/ODESystemSteppers.h"

#include "mml/algorithms/RootFinding.h"
#include "mml/algorithms/FunctionsAnalyzer.h"

#include "mpl/Mechanics/Pendulums.h"
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

	vvp.print();

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
	// ComparePendulumPeriods(sys, initCond, t1, t2, expectNumSteps, minSaveInterval);  // Defined in pendulum_comparison.cpp

	// Comparing number of adaptive steps needed, for given accuracy
	// CompareAdaptiveStepCounts(sys, initCond, t1, t2, minSaveInterval);  // Defined in pendulum_comparison.cpp

	ComparePhaseSpace(sys, initCond, t1, t2, minSaveInterval);
}
