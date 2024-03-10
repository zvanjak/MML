#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/Vector.h"

#include "core/Visualizer.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemSteppers.h"
#endif

#include "../test_data/diff_eq_systems_test_bed.h"

using namespace MML;

void Demo1_Lorenz_multi_func()
{
    TestBeds::LorenzSystemODE sys0(10.0, 28.0, 8.0/3.0);
     
    const double atol=1.0e-3, rtol=atol, h1=0.01, hmin=0.0, x1=0.0, x2=50.0;

    std::cout << "\n***********************************************************\n";
    std::cout << "**********         Demo Lorenz system                ********\n";

    Vector<Real> ystart0(sys0.getDim());
    ystart0[0]=2.0;
    ystart0[1]=1.0;
    ystart0[2]=1.0;
    Output out0(10000);
    
    ODESystemSolver<StepperDopr5> ode_solver0(sys0,atol,rtol, out0);
    ODESystemSolution             sol0 = ode_solver0.integrate(ystart0, x1, x2, h1, hmin);

    Visualizer::VisualizeODESysSolAsMultiFunc(sol0, "Lorenz system", "demo1_lorenz_system.txt");
}

