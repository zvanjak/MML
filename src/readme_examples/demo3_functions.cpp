#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Vector.h"
#include "core/Matrix.h"

#include "core/Function.h"
#endif

#include "../test_data/real_functions_test_bed.h"
#include "../test_data/scalar_functions_test_bed.h"
#include "../test_data/vector_functions_test_bed.h"
#include "../test_data/parametric_curves_test_bed.h"

using namespace MML;

double Readme_functions_TestFunc(double x) 
{ 
    return sin(x)*(1.0 + 0.5*x*x); 
}
void Readme_functions()
{
    // creating a function object from an already existing (standalone) function
    RealFunction f1(Readme_functions_TestFunc);

    // or creating a function object directly
    RealFunction f2{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };

    // creating directly different types of functions
    ScalarFunction<3>       funcScalar([](const VectorN<Real, 3> &x) { return x[0]; });
    VectorFunction<3>       funcVector([](const VectorN<Real, 3> &x) { return VectorN<Real, 3>{0, x[0] * x[1], 0}; });
    VectorFunctionNM<2, 3>  funcVectorNM([](const VectorN<Real, 2> &x) { return VectorN<Real, 3>{0, x[0] * x[1], 0}; });
    ParametricCurve<3>      paramCurve([](double x) { return VectorN<Real, 3>{x, 2 * x, 3 * x}; });
    ParametricSurface<3>    paramSurface([](double x, double y) { return VectorN<Real, 3>{x * y, 2 * x * y, 3 * x}; });   

    // using predefined functions from TestBeds
    auto fdef1 = TestBeds::RealFunctionsTestBed::getTestFunctionReal("Sin");
    auto fdef2 = TestBeds::ScalarFunctionsTestBed::getTestFunctionScalar3(0);
    auto fdef3 = TestBeds::VectorFunctionsTestBed::getTestFunctionVector(0);
    auto fdef4 = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix");
}
