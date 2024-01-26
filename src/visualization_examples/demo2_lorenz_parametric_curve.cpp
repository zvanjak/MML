#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/Vector.h"
#include "core/Curves.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemSteppers.h"
#endif

#include "../test_data/diff_eq_systems_test_bed.h"

using namespace MML;

void Demo1_Lorenz_parametric_curve()
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

    // TODO - izvuci kao parametric curve rjesenje diff sustava

    // sol0.Serialize("lorenz_1.txt", "Lorenz system");
    // auto ret2 = std::system("..\\..\\tools\\visualizers\\real_function_visualizer\\MML_RealFunctionVisualizer.exe lorenz_1.txt");

    Curves::HelixCurve              helix(20.0, 20.0);
    Curves::ToroidalSpiralCurve     toroid(20.0);

    // helix.SerializeCartesian3D(0.0, 2.0 * Constants::PI, 100, "helix.txt");
    // std::system("..\\..\\tools\\visualizers\\parametric_curve_visualizer\\MML_ParametricCurveVisualizer.exe helix.txt");

    toroid.SerializeCartesian3D(0.0, 5.0 * Constants::PI, 500, "..\\..\\results\\toroid.txt");
    auto ret = std::system("..\\..\\tools\\visualizers\\parametric_curve_visualizer\\MML_ParametricCurveVisualizer.exe ..\\..\\results\\toroid.txt");    
}

