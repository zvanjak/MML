#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "utilities/StdFunctions.h"

#include "core/VectorN.h"
#include "core/Function.h"

#include "basic_types/Curves.h"
#include "basic_types/Surfaces.h"
#endif

using namespace MML;

double Demo_Function_TestFunc(double x) 
{ 
    return sin(x)*(1.0 + 0.5*x*x); 
}

void Demo_Function_from_funct_ptr()
{
    // creating a function object from a already existing (standalone) function
    RealFunction f1(Demo_Function_TestFunc);

    // or creating a function object directly
    RealFunction f2{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };
}

// TODO - finish Demo_Function_from_std_func
void Demo_Function_from_std_func()
{

}

void Demo_Function()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                          FUNCTION                             ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    std::function<double(double)> f(Demo_Function_TestFunc);
    
    RealFunction f2{[](double x) { return sin(x); } };
    RealFunctionFromStdFunc func(f);

    RealFunction func2([](double x) { return x; });
    ScalarFunction<3> funcScalar([](const VectorN<Real, 3> &x) { return x[0]; });
    VectorFunction<3> funcVector([](const VectorN<Real, 3> &x) { return VectorN<Real, 3>{0, x[0] * x[1], 0}; });

    ParametricCurve<3> paramCurve([](double x) { return VectorN<Real, 3>{x, 2 * x, 3 * x}; });

    //auto val1 = Curves::helix_curve(1.0);
    auto val2 = Surfaces::test1(1.0, 1.0);

    f2.SerializeVariableSpaced(0.0, 10.0, 100, "test.txt");
    auto ret = std::system("..\\..\\tools\\visualizers\\real_function_visualizer\\MML_RealFunctionVisualizer.exe");
}