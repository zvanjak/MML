#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Visualizer.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemSteppers.h"
#endif

#include "../test_data/diff_eq_systems_defs.h"
#include "../test_data/diff_eq_systems_test_bed.h"

using namespace MML;


void Readme_ode_solvers()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                    README - ODE solvers                       ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// in-place definition of Lorenz system (have to fix parameters!)
	ODESystem lorenzSystem(3, [](Real t, const Vector<Real>& y, Vector<Real>& dydt)
		{
			double sigma = 10.0, rho = 28.0, beta = 8.0 / 3.0;
			dydt[0] = sigma * (y[1] - y[0]);
			dydt[1] = y[0] * (rho - y[2]) - y[1];
			dydt[2] = y[0] * y[1] - beta * y[2];
		});

	// get it from predefined test-bed 
	TestBeds::LorenzSystemODE alsoLorenzSystem(10.0, 28.0, 8.0 / 3.0);

	const double atol = 1.0e-3, rtol = atol, h1 = 0.01, hmin = 0.0;
	double x1 = 0.0, x2 = 50.0;

	Vector<Real> init_cond({ 2.0, 1.0, 1.0 });
	Output out0(10000);

	ODESystemSolver<StepperDopr5> ode_solver0(lorenzSystem, atol, rtol, out0);
	ODESystemSolution             sol0 = ode_solver0.integrate(init_cond, x1, x2, h1, hmin);

	//Visualizer::VisualizeODESysSolAsMultiFunc(sol0, "Lorenz system", "demo5_lorenz_system.txt");
	Visualizer::VisualizeODESysSolAsParamCurve3(sol0, "Lorenz system as parametric curve", "demo5_lorenz_system_as_parametric_curve.txt");

	Vector<Real> init_cond1({ 25.0, 15.0, -15.0 });
	Vector<Real> init_cond2({ -20.0, -10.0, -5.0 });
	Vector<Real> init_cond3({ -25.0, 50.0, -10.0 });
	Vector<Real> init_cond4({ 20.0, -20.0, 30.0 });
	Vector<Real> init_cond5({ -70.0, 20.0, -50.0 });

	Output out1(10000);
	ODESystemSolver<StepperDopr5> ode_solver1(lorenzSystem, atol, rtol, out1);
	ODESystemSolution             sol1 = ode_solver1.integrate(init_cond1, x1, x2, h1, hmin);
	Visualizer::VisualizeODESysSolAsParamCurve3(sol1, "Lorenz system as parametric curve", "demo5_lorenz_system_as_parametric_curve_1.txt");

	Output out2(10000);
	ODESystemSolver<StepperDopr5> ode_solver2(lorenzSystem, atol, rtol, out2);
	ODESystemSolution             sol2 = ode_solver2.integrate(init_cond2, x1, x2, h1, hmin);
	Visualizer::VisualizeODESysSolAsParamCurve3(sol2, "Lorenz system as parametric curve", "demo5_lorenz_system_as_parametric_curve_2.txt");

	Output out3(10000);
	ODESystemSolver<StepperDopr5> ode_solver3(lorenzSystem, atol, rtol, out3);
	ODESystemSolution             sol3 = ode_solver3.integrate(init_cond3, x1, x2, h1, hmin);
	Visualizer::VisualizeODESysSolAsParamCurve3(sol3, "Lorenz system as parametric curve", "demo5_lorenz_system_as_parametric_curve_3.txt");

	Output out4(10000);
	ODESystemSolver<StepperDopr5> ode_solver4(lorenzSystem, atol, rtol, out4);
	ODESystemSolution             sol4 = ode_solver4.integrate(init_cond4, x1, x2, h1, hmin);
	Visualizer::VisualizeODESysSolAsParamCurve3(sol4, "Lorenz system as parametric curve", "demo5_lorenz_system_as_parametric_curve_4.txt");

	Output out5(10000);
	ODESystemSolver<StepperDopr5> ode_solver5(lorenzSystem, atol, rtol, out5);
	ODESystemSolution             sol5 = ode_solver5.integrate(init_cond5, x1, x2, h1, hmin);
	Visualizer::VisualizeODESysSolAsParamCurve3(sol5, "Lorenz system as parametric curve", "demo5_lorenz_system_as_parametric_curve_5.txt");

	// visualize all together
	Visualizer::VisualizeMultiParamCurve3D({ 
										"demo5_lorenz_system_as_parametric_curve_1.txt",
										"demo5_lorenz_system_as_parametric_curve_2.txt",
										"demo5_lorenz_system_as_parametric_curve_3.txt",
										"demo5_lorenz_system_as_parametric_curve_4.txt",
										"demo5_lorenz_system_as_parametric_curve_5.txt"
									});
}