#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Visualizer.h"
#include "core/Serializer.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemSteppers.h"
#include "algorithms/ODESystemSteppers_Stiff.h"
#endif

#include "../test_data/diff_eq_systems_test_bed.h"

using namespace MML;

void Docs_Demo_ODE_solvers_5th_order()
{
	std::cout << "\n***********************************************************\n";
	std::cout << "***         Solving Lorenz system - 5th order             ***\n";

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

	//std::cout << "x values:\n";
	//sol0._xval.Print(std::cout, 6, 3); std::cout << std::endl;
	//std::cout << "y values: - ";
	//sol0._yval.Print(std::cout, 6, 3);

	Visualizer::VisualizeODESysSolAsMultiFunc(sol0, "Lorenz system as multi function", "lorenz_1_as_multi_func.txt");
	Visualizer::VisualizeODESysSolAsParamCurve3(sol0, "Lorenz system as parametric curve", "demo5_lorenz_system_as_parametric_curve.txt");
}

void Docs_Demo_ODE_solvers_8th_order()
{
	std::cout << "\n***********************************************************\n";
	std::cout << "***         Solving Lorenz system - 8th order             ***\n";
	// get it from predefined test-bed 
	TestBeds::LorenzSystemODE sys0(10.0, 28.0, 8.0 / 3.0);

	const double atol = 1.0e-3, rtol = atol, h1 = 0.01, hmin = 0.0;
	double x1 = 0.0, x2 = 50.0;
	Vector<Real> ystart01({ 2.0, 1.0, 1.0 });
	Output out01(10000);

	ODESystemSolver<StepperDopr853> ode_solver01(sys0, atol, rtol, out01);
	ODESystemSolution               sol01 = ode_solver01.integrate(ystart01, x1, x2, h1, hmin);

	Visualizer::VisualizeODESysSolAsMultiFunc(sol01, "Lorenz system solution - Dormand 8th order", "lorenz_dormand8.txt");
}

void Docs_Demo_ODE_solvers_BulirschStoer()
{
	std::cout << "\n***********************************************************\n";
	std::cout << "**********            Bulirsch-Stoer               ********\n";
	// get it from predefined test-bed 
	TestBeds::LorenzSystemODE sys0(10.0, 28.0, 8.0 / 3.0);

	const double atol = 1.0e-3, rtol = atol, h1 = 0.01, hmin = 0.0;
	double x1 = 0.0, x2 = 50.0;
	Vector<Real> ystartBS({ 2.0, 1.0, 1.0 });
	Output out_BS(10000);

	ODESystemSolver<StepperBS> ode_solver_BS(sys0, atol, rtol, out_BS);
	ODESystemSolution          sol_BS = ode_solver_BS.integrate(ystartBS, x1, x2, h1, hmin);

	//std::cout << "x values:\n";    sol_BS._xval.Print(std::cout, 6, 3); std::cout << std::endl;
	//std::cout << "y values: - ";   sol_BS._yval.Print(std::cout, 6, 3);

	Visualizer::VisualizeODESysSolAsMultiFunc(sol_BS, "Lorenz system solution - Bulirsch-Stoer", "lorenz_bulirsch_stoer.txt");
}

void Docs_Demo_ODE_solvers_RungeKutta_4th_order()
{
	std::cout << "\n***********************************************************\n";
	std::cout << "******        Runge-Kutta 4th order - Dumb           ******\n";
	// get it from predefined test-bed 
	TestBeds::LorenzSystemODE sys0(10.0, 28.0, 8.0 / 3.0);

	const double atol = 1.0e-3, rtol = atol, h1 = 0.01, hmin = 0.0;
	double x1 = 0.0, x2 = 50.0;
	Vector<Real> ystart0({ 2.0, 1.0, 1.0 });
	//Output out_BS(10000);

	RungeKuttaSolverDumb    rkdumb;
	ODESystemSolutionEqualSpacing sol = rkdumb.integrate(sys0, ystart0, 0.0, 2.0, 200);

	std::cout << "x values:\n";    sol.xval.Print(std::cout, 6, 3); std::cout << std::endl;
	std::cout << "y values: - ";   sol.yval.Print(std::cout, 6, 3);

	std::cout << "\n***********************************************************\n";
	std::cout << "******       Runge-Kutta 4th order - stepper         ******\n";

	int nok, nbad;
	ystart0[0] = 2.0;
	ystart0[1] = 1.0;
	ystart0[2] = 1.0;
	RungeKuttaSolverSimple   rkNR2;
	ODESystemSolution sol2 = rkNR2.integrate(sys0, ystart0, 0.0, 2.0, 100, 0.1, 1e-06, h1, hmin, nok, nbad);

	std::cout << "x values:\n";    sol2._xval.Print(std::cout, 7, 3); std::cout << std::endl;
	std::cout << "y values: - ";   sol2._yval.Print(std::cout, 7, 3);

	Visualizer::VisualizeODESysSolAsMultiFunc(sol2, "Lorenz system solution - RK 4th order", "lorenz_RK_4th.txt");
}

void Docs_Demo_ODE_solvers_Stiff_Rosenbrock()
{
	std::cout << "\n***********************************************************\n";
	std::cout << "******      Stiff system - Rosenbrock method         ******\n";

	auto sys_stiff = TestBeds::ODESystemTestBed::getODESystemWithJacobian(0);

	const double atol = 1.0e-3, rtol = atol, h1 = 0.01, hmin = 0.0, x1 = 0.0, x2 = 2.0;
	Vector<Real> ystart02(sys_stiff.getDim());
	ystart02[0] = 1.0;
	ystart02[1] = 1.0;
	ystart02[2] = 0.0;
	Output out02(20);

	ODESystemSolver<StepperRoss> ode_solver02(sys_stiff, atol, rtol, out02);
	ODESystemSolution            sol02 = ode_solver02.integrate(ystart02, 0.0, 50.0, h1, hmin);

	std::cout << "x values:\n";    sol02._xval.Print(std::cout, 6, 3); std::cout << std::endl;
	std::cout << "y values: - ";   sol02._yval.Print(std::cout, 6, 3);
}

void Docs_Demo_ODE_solvers_Stiff_SemiImplExtrapolation_method()
{
	std::cout << "\n***********************************************************\n";
	std::cout << "****    Stiff system - Semi-implicit extrapol method    ***\n";

	auto sys_stiff = TestBeds::ODESystemTestBed::getODESystemWithJacobian(0);

	const double atol = 1.0e-3, rtol = atol, h1 = 0.01, hmin = 0.0, x1 = 0.0, x2 = 2.0;
	Vector<Real> ystart03(sys_stiff.getDim());
	ystart03[0] = 1.0;
	ystart03[1] = 1.0;
	ystart03[2] = 0.0;
	Output out03(20);

	ODESystemSolver<StepperSemiImplExtr> ode_solver03(sys_stiff, atol, rtol, out03);
	ODESystemSolution        sol03 = ode_solver03.integrate(ystart03, 0.0, 50.0, h1, hmin);

	std::cout << "x values:\n";    sol03._xval.Print(std::cout, 6, 3); std::cout << std::endl;
	std::cout << "y values: - ";   sol03._yval.Print(std::cout, 6, 3);
}

void Docs_Demo_ODE_solvers()
{
	Docs_Demo_ODE_solvers_5th_order();
	Docs_Demo_ODE_solvers_8th_order();
	Docs_Demo_ODE_solvers_BulirschStoer();
	Docs_Demo_ODE_solvers_RungeKutta_4th_order();
	Docs_Demo_ODE_solvers_Stiff_Rosenbrock();
	Docs_Demo_ODE_solvers_Stiff_SemiImplExtrapolation_method();
}