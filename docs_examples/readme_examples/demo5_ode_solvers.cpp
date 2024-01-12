#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemSteppers.h"
#endif

#include "../test_data/diff_eq_systems_test_bed.h"

using namespace MML;


void Readme_ode_solvers()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                    README - ODE solvers                       ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // in-place definition
    ODESystem vanDerPol(2, [](double t, const Vector<Real>& y, Vector<Real>& dydt)
    {
        double mju = 0.1;
        dydt[0] = y[1];
        dydt[1] = mju * (1.0 - y[0] * y[0]) * y[1] - y[0];
    });
    // get it from predefined test-bed (but, it has a fixed parameter value!)
    auto vanDerPolAlso  = TestBeds::ODESystemTestBed::getODESystem("VanDerPol 0.1");

    const double atol=1.0e-3, rtol=atol;
    const double h1=0.01, hmin=0.0;
    const double x1=0.0, x2=2.0;
    
    Vector<Real> init_cond{2.0, 0.0};
    Output out(20);      // saving only 20 points     

    ODESystemSolver<StepperDopr853> ode_solver(vanDerPol, atol, rtol, out);
    ODESystemSolution               sol = ode_solver.integrate(init_cond, x1, x2, h1, hmin);

    std::cout << "x values:\n";  sol._xval.Print(std::cout, 6, 3); std::cout << std::endl;
    std::cout << "y values: - "; sol._yval.Print(std::cout, 6, 3);
/* OUTPUT
x values:
[     0,    0.1,    0.2,    0.3,    0.4,    0.5,    0.6,    0.7,    0.8,    0.9,      1,    1.1,    1.2,    1.3,    1.4,    1.5,    1.6,    1.7,    1.8,    1.9,      2]
y values: - Rows: 2 Cols: 21
[      2,   1.99,   1.96,   1.91,   1.85,   1.77,   1.67,   1.56,   1.43,   1.29,   1.14,  0.976,  0.804,  0.623,  0.436,  0.242, 0.0439, -0.157, -0.358, -0.557, -0.752,  ]
[      0, -0.197, -0.386, -0.567, -0.738, -0.901,  -1.05,   -1.2,  -1.33,  -1.45,  -1.57,  -1.67,  -1.77,  -1.85,  -1.91,  -1.96,     -2,  -2.01,     -2,  -1.97,  -1.92,  ]
*/
}