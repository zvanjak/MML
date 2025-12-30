#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/Vector.h"

#include "tools/Visualizer.h"
#include "tools/Serializer.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/ODESystemSteppers.h"
#endif

#include "../test_data/diff_eq_systems_test_bed.h"

using namespace MML;

void Demo_Lorenz_solve()
{
	//TestBeds::LorenzSystemODE sys0(10.0, 28.0, 8.0/3.0);
	// 
	//const double atol=1.0e-3, rtol=atol, h1=0.01, hmin=0.0, x1=0.0, x2=50.0;

	//std::cout << "\n***********************************************************\n";
	//std::cout << "**********         Demo Lorenz system                ********\n";

	//Vector<Real> ystart0({2.0, 1.0, 1.0});
	//Output out0(10000);
	//
	//ODESystemSolver<StepperDopr5> ode_solver0(sys0,atol,rtol, out0);
	//ODESystemSolution             sol0 = ode_solver0.integrate(ystart0, x1, x2, h1, hmin);

	//std::cout << "x values:\n";
	//sol0._tval.Print(std::cout, 6, 3); std::cout << std::endl;
	//std::cout << "y values: - ";
	//sol0._xval.Print(std::cout, 6, 3);

	//Visualizer::VisualizeODESysSolAsMultiFunc(sol0, "Lorenz system as multi function", "lorenz_1_as_multi_func.txt");
}

void Demo_VanderPol_solve()
{
	auto sys0 = TestBeds::ODESystemTestBed::getODESystem(1);

	const double atol = 1.0e-3, rtol = atol, h1 = 0.01, hmin = 0.0, x1 = 0.0, x2 = 2.0;

	std::cout << "\n***********************************************************\n";
	std::cout << "******        Runge-Kutta 4th order - Dumb           ******\n";

	Vector<Real> ystart0(sys0.getDim());
	ystart0[0] = 2.0;
	ystart0[1] = 0.0;

	RungeKutta4_StepCalculator  rk4FixedStepCalc;
	ODESystemFixedStepSolver	solver1(sys0, rk4FixedStepCalc);
	ODESystemSolution sol = solver1.integrate(ystart0, 0.0, 2.0, 20);

	std::cout << "x values:\n";    sol.getTValues().Print(std::cout, 6, 3); std::cout << std::endl;
	std::cout << "y values: - ";   sol.getXValues().Print(std::cout, 6, 3);

	std::cout << "\n***********************************************************\n";
	std::cout << "******       Runge-Kutta 4th order - stepper         ******\n";

	ystart0[0] = 2.0;
	ystart0[1] = 0.0;
	ODESystemSolver<RK5_CashKarp_Stepper> solver2(sys0);
	ODESystemSolution sol2 = solver2.integrate(ystart0, 0.0, 2.0, 0.1, 1e-06, h1, hmin);

	std::cout << "x values:\n";    sol2.getTValues().Print(std::cout, 7, 3); std::cout << std::endl;
	std::cout << "y values: - ";   sol2.getXValues().Print(std::cout, 7, 3);

	Visualizer::VisualizeODESysSolAsMultiFunc(sol2, "Van der Pol solution - RK 4th order", std::vector<std::string>{"comp 1", "comp 2"}, "vanderpol_RK_4th.txt");

}

void Demo_SimpleLinearODE_solve()
{
	//const double atol = 1.0e-3, rtol = atol, h1 = 0.01, hmin = 0.0, x1 = 0.0, x2 = 2.0;

	//auto sys = TestBeds::ODESystemTestBed::getTestODESystemWithSolution(0);

	//RungeKutta4_StepCalculator   rk4FixedStepCalc;
	//ODESystemFixedStepSolver  solver(*sys.getODESystem(), rk4FixedStepCalc);
	//ODESystemSolution sol = solver.integrate(sys.getInitialConditions(), 0.0, 2.0, 20);

	//std::cout << "**********    Runge-Kutta 4th order - Dumb    ********\n";
	//std::cout << "x values:\n";    sol.getTValues().Print(std::cout, 7, 3); std::cout << std::endl;
	//std::cout << "y values: - ";   sol.getXValues().Print(std::cout, 7, 3);

	//ODESystemSolver<RK5_CashKarp_Stepper> solver2(*sys.getODESystem());
	//ODESystemSolution sol2 = solver2.integrate(sys.getInitialConditions(), 0.0, 2.0, 0.1, 1e-06, h1, hmin);

	//std::cout << "**********    Runge-Kutta 4th order - stepper    ********\n";
	//std::cout << "x values:\n";    sol2.getTValues().Print(std::cout, 7, 3); std::cout << std::endl;
	//std::cout << "y values: - ";   sol2.getXValues().Print(std::cout, 7, 3);
}


void Demo_ODESystemSolvers()
{
	std::cout << std::endl;
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                    DIFF.EQUATIONS SOLVERS                     ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// Demo_VanderPol_solve();

	// Demo_Lorenz_solve();

	// int _numPoints = TestBeds::ODESystemTestBed::numODESystem();
	// Demo_SimpleLinearODE_solve();
}