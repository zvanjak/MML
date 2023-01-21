#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include <iostream>
#include <iomanip>
#include <cmath>

#include "basic_types/VectorN.h"
#include "basic_types/Function.h"
#endif

double TestFunc1(double x)
{
    return x;
}

void Demo_Function()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                          FUNCTION                             ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    std::function<double(double)> f(TestFunc1);
    
    MML::RealFunctionFromStdFunc func(f);

    MML::RealFunction func2([](double x) { return x; });
    MML::ScalarFunctionFromFuncPtr<3> funcScalar([](const MML::VectorN<Real, 3> &x) { return x[0]; });
    MML::VectorFunctionFromFuncPtr<3> funcVector([](const MML::VectorN<Real, 3> &x) { return MML::VectorN<Real, 3>{0, x[0] * x[1], 0}; });

    MML::ParametricCurveFromFuncPtr<3> paramCurve([](double x) { return MML::VectorN<Real, 3>{x, 2 * x, 3 * x}; });
}