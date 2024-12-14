#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "tools/Visualizer.h"
#include "tools/Serializer.h"

#include "algorithms/ODESystemSolver.h"
#endif

#include "../test_data/diff_eq_systems_test_bed.h"

using namespace MML;

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


void Docs_Demo_ODE_solvers()
{
	Docs_Demo_ODE_solvers_RungeKutta_4th_order();
}