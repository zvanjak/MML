#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Vector.h"
#include "core/Matrix.h"

#include "core/Function.h"
#include "basic_types/Curves.h"
#include "basic_types/Geometry3D.h"

#include "algorithms/ODESystemSolvers.h"
#include "algorithms/ODESystemSteppers.h"
#include "algorithms/ODESystemSteppers_Stiff.h"
#include "algorithms/DiffGeometryAlgorithms.h"
#endif

#include "../test_data/parametric_curves_test_bed.h"
#include "../test_data/diff_eq_systems_test_bed.h"
#include "../test_data/linear_alg_eq_systems_test_bed.h"

using namespace MML;


void Readme_ode_solvers()
{
    auto sys0 = TestBeds::ODESystemTestBed::getODESystem(1);

    const double atol=1.0e-3, rtol=atol;
    const double h1=0.01, hmin=0.0;
    const double x1=0.0, x2=2.0;
    
    Vector<Real> ystart01{2.0, 0.0};
    Output out(20);            

    ODESystemSolver<StepperDopr853> ode_solver01(sys0, atol, rtol, out);
    ODESystemSolution               sol01 = ode_solver01.integrate(ystart01, 0.0, 2.0, h1, hmin);

    std::cout << "x values:\n";  sol01._xval.Print(std::cout, 6, 3); std::cout << std::endl;
    std::cout << "y values: - "; sol01._yval.Print(std::cout, 6, 3);
}