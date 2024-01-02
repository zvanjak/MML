#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Vector.h"
#include "core/Matrix.h"

#include "core/Function.h"
#include "core/Derivation.h"
#include "core/Integration.h"
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

    // numerical derivation of real function (different orders)
    double der_f1 = Derivation::NDer1(f2, 0.5);
    double der_f2 = Derivation::NDer2(f1, 0.5);
    double der_f4 = Derivation::NDer4(f1, 0.5, 1e-6);   // setting explicit step size
    double der_f6 = Derivation::NDer6(f2, 0.5);    
    double der_f8 = Derivation::NDer8(f2, 0.5);
    // we can use default Derive routine (default set to NDer4)
    double num_der4 = Derivation::Derive(f1, 0.5, nullptr);

    double a = 0.0;
    double b = 1.0;
    double int_trap = Integration::IntegrateTrap(f1,a,b);
    double int_simp = Integration::IntegrateSimpson(f1,a,b);
    double int_romb = Integration::IntegrateRomberg(f1,a,b);
    // we can use default Integrate routine (default set to IntegrateSimpson)
    double int_def = Integration::Integrate(f1, a, b, 1e-04);        
}
