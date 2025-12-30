#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "mml/MMLBase.h"

#include "mml/base/Vector.h"

#include "mml/tools/Visualizer.h"

#include "mml/algorithms/ODESystemSolver.h"
#include "mml/algorithms/ODESystemStepCalculators.h"

#include "mpl/Mechanics/HarmonicOscillator.h"
#include "mpl/Mechanics/Oscillators.h"
#endif

using namespace MML;
using namespace MPL;

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
