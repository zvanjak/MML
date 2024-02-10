#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

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
    ODESystem LorenzSystem(3, [](Real t, const Vector<Real>& y, Vector<Real>& dydt)
    {
        double sigma = 10.0, rho = 28.0, beta = 8.0 / 3.0;
        dydt[0] = sigma * (y[1] - y[0]);
        dydt[1] = y[0] * (rho - y[2]) - y[1];
        dydt[2] = y[0] * y[1] - beta * y[2];
    });    

    // get it from predefined test-bed 
    TestBeds::LorenzSystemODE alsoLorenzSystem(10.0, 28.0, 8.0/3.0);

    const double atol=1.0e-3, rtol=atol, h1=0.01, hmin=0.0;
    double x1=0.0, x2=50.0;
    
    Vector<Real> init_cond({2.0, 1.0, 1.0});
    Output out0(10000);
    
    ODESystemSolver<StepperDopr5> ode_solver0(LorenzSystem,atol,rtol, out0);
    ODESystemSolution             sol0 = ode_solver0.integrate(init_cond, x1, x2, h1, hmin);

    sol0.Serialize("..\\..\\results\\demo5_lorenz_system.txt", "Lorenz system");
    std::system("..\\..\\tools\\visualizers\\real_function_visualizer\\MML_RealFunctionVisualizer.exe ..\\..\\results\\demo5_lorenz_system.txt");

    auto curve = sol0.getSolutionAsParametricCurve<3>();
    sol0.SerializeAsParametricCurve3D("..\\..\\results\\demo5_lorenz_system_as_parametric_curve.txt", "Lorenz system as parametric curve");
    std::system("..\\..\\tools\\visualizers\\parametric_curve_visualizer\\MML_ParametricCurveVisualizer.exe ..\\..\\results\\demo5_lorenz_system_as_parametric_curve.txt");
}