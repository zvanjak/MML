#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "utilities/StdFunctions.h"

#include "base/VectorN.h"
#include "core/Function.h"

#include "core/Curves.h"
#include "core/Surfaces.h"
#endif

using namespace MML;

Real Demo_Function_TestFunc(Real x) 
{ 
    return sin(x)*(1 + x*x / 2); 
}

void Demo_Function_from_funct_ptr()
{
    // creating a function object from a already existing (standalone) function
    RealFunction f1(Demo_Function_TestFunc);

    // or creating a function object directly
    RealFunction f2{[](Real x) { return sin(x)*(1 + x*x / 2); } };
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

    std::function<Real(Real)> f(Demo_Function_TestFunc);
    
    RealFunction f1(Demo_Function_TestFunc);
    RealFunction f2{ [](Real x) { return sin(x); } };
    RealFunctionFromStdFunc f3(f);

    ScalarFunction<3>       funcScalar([](const VectorN<Real, 3> &x) { return x[0]; });
    VectorFunction<3>       funcVector([](const VectorN<Real, 3> &x) { return VectorN<Real, 3>{0, x[0] * x[1], 0}; });
    VectorFunctionNM<2, 3>  funcVectorNM([](const VectorN<Real, 2> &x) { return VectorN<Real, 3>{0, x[0] * x[1], 0}; });
    ParametricCurve<3>      paramCurve([](Real x) { return VectorN<Real, 3>{x, 2 * x, 3 * x}; });
    ParametricSurface<3>    paramSurface([](Real x, Real y) { return VectorN<Real, 3>{x * y, 2 * x * y, 3 * x}; });

    auto val2 = Surfaces::test1(1.0, 1.0);

    VectorN<Real, 3> p1{2.0, 2.0, 5};
    
    MatrixNM<Real, 3, 3> jac = funcVector.jacobian(p1);    
}