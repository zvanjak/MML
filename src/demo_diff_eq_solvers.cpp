#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Vector.h"

#include "algorithms/ODESystemSolvers.h"
#include "algorithms/ODESystemSteppers.h"
#include "algorithms/ODESystemSteppers_Stiff.h"
#endif

#include "../test_data/diff_eq_systems_test_bed.h"

using namespace MML;


void Demo_VanderPol_solve()
{
    auto sys0 = TestBeds::ODESystemTestBed::getODESystem(1);
    auto sys_stiff = TestBeds::ODESystemTestBed::getODESystemWithJacobian(0);

    const double atol=1.0e-3, rtol=atol, h1=0.01, hmin=0.0, x1=0.0, x2=2.0;

    std::cout << "\n***********************************************************\n";
    std::cout << "**********         Dopler 5th order                ********\n";

    Vector<Real> ystart0(sys0.getDim());
    ystart0[0]=2.0;
    ystart0[1]=0.0;
    Output out0(20);         
    
    ODESystemSolver<StepperDopr5> ode_solver0(sys0,atol,rtol, out0);
    ODESystemSolution             sol0 = ode_solver0.integrate(ystart0, x1, x2, h1, hmin);

    std::cout << "x values:\n";
    sol0.xval.Print(std::cout, 6, 3); std::cout << std::endl;
    std::cout << "y values: - ";
    sol0.yval.Print(std::cout, 6, 3);

    std::cout << "\n***********************************************************\n";
    std::cout << "**********         Dopler 8th order                ********\n";
    
    Vector<Real> ystart01(sys0.getDim());
    ystart01[0]=2.0;
    ystart01[1]=0.0;
    Output out01(20);            

    ODESystemSolver<StepperDopr853> ode_solver01(sys0, atol,rtol,out01);
    ODESystemSolution             sol01 = ode_solver01.integrate(ystart01, 0.0, 2.0, h1, hmin);

    std::cout << "x values:\n";    sol01.xval.Print(std::cout, 6, 3); std::cout << std::endl;
    std::cout << "y values: - ";   sol01.yval.Print(std::cout, 6, 3);

    std::cout << "\n***********************************************************\n";
    std::cout << "**********            Bulirsch-Stoer               ********\n";
    
    Vector<Real> ystartBS(sys0.getDim());
    ystartBS[0]=2.0;
    ystartBS[1]=0.0;
    Output out_BS(20);            

    ODESystemSolver<StepperBS> ode_solver_BS(sys0, atol,rtol,out_BS);
    ODESystemSolution          sol_BS = ode_solver_BS.integrate(ystartBS, 0.0, 2.0, h1, hmin);

    std::cout << "x values:\n";    sol_BS.xval.Print(std::cout, 6, 3); std::cout << std::endl;
    std::cout << "y values: - ";   sol_BS.yval.Print(std::cout, 6, 3);
 
    std::cout << "\n***********************************************************\n";
    std::cout << "******        Runge-Kutta 4th order - Dumb           ******\n";    
    
    ystart0[0]=2.0;
    ystart0[1]=0.0;
    RungeKuttaSolverDumb    rkdumb;
    ODESystemSolutionEqualSpacing sol = rkdumb.integrate(sys0, ystart0, 0.0, 2.0, 20);

    std::cout << "x values:\n";    sol.xval.Print(std::cout, 6, 3); std::cout << std::endl;
    std::cout << "y values: - ";   sol.yval.Print(std::cout, 6, 3);

    std::cout << "\n***********************************************************\n";
    std::cout << "******       Runge-Kutta 4th order - stepper         ******\n";  

    int nok,  nbad;
    ystart0[0]=2.0;
    ystart0[1]=0.0;
    RungeKuttaNR2    rkNR2;
    ODESystemSolution sol2 = rkNR2.integrate(sys0, ystart0, 0.0, 2.0, 20, 0.1, 1e-06, h1, hmin, nok, nbad);

    std::cout << "x values:\n";    sol2.xval.Print(std::cout, 7,3); std::cout << std::endl;
    std::cout << "y values: - ";   sol2.yval.Print(std::cout, 7,3); 


    std::cout << "\n***********************************************************\n";
    std::cout << "******      Stiff system - Rosenbrock method         ******\n";

    Vector<Real> ystart02(sys_stiff.getDim());
    ystart02[0]=1.0;
    ystart02[1]=1.0;
    ystart02[2]=0.0;
    Output out02(20);            

    ODESystemSolver<StepperRoss> ode_solver02(sys_stiff, atol, rtol, out02);
    ODESystemSolution            sol02 = ode_solver02.integrate(ystart02, 0.0, 50.0, h1, hmin);

    std::cout << "x values:\n";    sol02.xval.Print(std::cout, 6, 3); std::cout << std::endl;
    std::cout << "y values: - ";   sol02.yval.Print(std::cout, 6, 3);
    
    std::cout << "\n***********************************************************\n";
    std::cout << "****    Stiff system - Semi-implicit extrapol method    ***\n";

    Vector<Real> ystart03(sys_stiff.getDim());
    ystart03[0]=1.0;
    ystart03[1]=1.0;
    ystart03[2]=0.0;
    Output out03(20);            

    ODESystemSolver<StepperSemiImplExtr> ode_solver03(sys_stiff, atol, rtol, out03);
    ODESystemSolution        sol03 = ode_solver03.integrate(ystart03, 0.0, 50.0, h1, hmin);

    std::cout << "x values:\n";    sol03.xval.Print(std::cout, 6, 3); std::cout << std::endl;
    std::cout << "y values: - ";   sol03.yval.Print(std::cout, 6, 3);              
}

void Demo_SimpleLinearODE_solve()
{
    const double atol=1.0e-3, rtol=atol, h1=0.01, hmin=0.0, x1=0.0, x2=2.0;

    auto sys = TestBeds::ODESystemTestBed::getTestODESystemWithSolution(0);

    RungeKuttaSolverDumb    rkdumb;
    ODESystemSolutionEqualSpacing sol = rkdumb.integrate(sys, sys.getInitialConditions(), 0.0, 2.0, 20);

    std::cout << "**********    Runge-Kutta 4th order - Dumb    ********\n";
    std::cout << "x values:\n";    sol.xval.Print(std::cout, 7, 3); std::cout << std::endl;
    std::cout << "y values: - ";   sol.yval.Print(std::cout, 7, 3);
    
    RungeKuttaNR2    rkNR2;
    int nok,  nbad;
    ODESystemSolution sol2 = rkNR2.integrate(sys, sys.getInitialConditions(), 0.0, 2.0, 20, 0.1, 1e-06, h1, hmin, nok, nbad);

    std::cout << "**********    Runge-Kutta 4th order - stepper    ********\n";
    std::cout << "x values:\n";    sol2.xval.Print(std::cout, 7,3); std::cout << std::endl;
    std::cout << "y values: - ";   sol2.yval.Print(std::cout, 7,3);
}

void Demo_DiffEqSolvers()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                    DIFF.EQUATIONS SOLVERS                     ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Demo_VanderPol_solve();

    int n = TestBeds::ODESystemTestBed::numODESystem();
    Demo_SimpleLinearODE_solve();
}