#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Vector.h"
#include "core/Matrix.h"

#include "core/Function.h"
#include "core/InterpolatedFunction.h"
#endif

#include "../test_data/real_functions_test_bed.h"
#include "../test_data/scalar_functions_test_bed.h"
#include "../test_data/vector_functions_test_bed.h"
#include "../test_data/parametric_curves_test_bed.h"

using namespace MML;

/////////////////////////////////////////////////////////////////////////////////////
// CASE 1 - standalone function providing calculation of a function
double Readme_functions_TestFunc(double x) { 
    return sin(x)*(1.0 + 0.5*x*x); 
}

void Readme_defining_functions_case_1_usage()
{
    // creating a function object from an already existing (standalone) function
    RealFunction f1(Readme_functions_TestFunc);
}

/////////////////////////////////////////////////////////////////////////////////////
// CASE 2 - create it directly with lambda
void Readme_defining_functions_case_2_usage()
{
    RealFunction f2{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };

    // creating directly different types of functions
    ScalarFunction<3>       fScalar([](const VectorN<Real, 3> &x) { return x[0]; });
    VectorFunction<3>       fVector([](const VectorN<Real, 3> &x) { return VectorN<Real, 3>{0, x[0] * x[1], 0}; });
    VectorFunctionNM<2, 3>  fVectorNM([](const VectorN<Real, 2> &x) { return VectorN<Real, 3>{0, x[0] * x[1], 0}; });
    ParametricCurve<3>      paramCurve([](double x) { return VectorN<Real, 3>{x, 2 * x, 3 * x}; });
    ParametricSurface<3>    paramSurface([](double x, double y) { return VectorN<Real, 3>{x * y, 2 * x * y, 3 * x}; });  

    // using predefined functions from TestBeds
    auto fdef1 = TestBeds::RealFunctionsTestBed::getTestFunctionReal("Sin");
    auto fdef2 = TestBeds::ScalarFunctionsTestBed::getTestFunctionScalar3(0);
    auto fdef3 = TestBeds::VectorFunctionsTestBed::getTestFunctionVector(0);
    auto fdef4 = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix"); 
}

/////////////////////////////////////////////////////////////////////////////////////
// CASE 3 - class you CAN change has member function that does the calculation
// Option 1 - define operator() for your class and create RealFunctionFromStdFunc
class ClassProvidingFuncToDerive {
    public:
        double operator()(double x ) const { 
            return 1.0;     /* calculation using member variables */ 
        }
};
// Option 2 - inherit class from IRealFunction interface, and use the object itself as RealFunction
class ClassProvidingFuncToDerive2 : public IRealFunction {
    public:
        double operator()(double x ) const { 
            return 1.0;     /* calculation using member variables */ 
        }
};
void Readme_defining_functions_case_3_usage()
{
    ClassProvidingFuncToDerive   obj1;
    RealFunctionFromStdFunc f1(std::function<double(double)>{obj1});
    
    ClassProvidingFuncToDerive2   f2;       // usable RealFunction object
}

/////////////////////////////////////////////////////////////////////////////////////
// CASE 4 - class you CAN'T change has member function that does the calculation
class BigComplexClassYouCantChange {
    // has data for calculating function you want to do somethign with
};

// Create a helper wrapper class, inherit it from IRealFunction and use it as RealFunction
class BigComplexFunc2 : public IRealFunction {
    const BigComplexClassYouCantChange &_ref;
public:
    BigComplexFunc2(const BigComplexClassYouCantChange &bigClass) : _ref(bigClass) { }

    double operator()(double x ) const {
        return 1.0;     /* calculation using _ref */ 
    }
};
void Readme_defining_functions_case_4_usage()
{
    BigComplexClassYouCantChange bigObj;

    BigComplexFunc2    f1(bigObj);     // usable RealFunction object
}

/////////////////////////////////////////////////////////////////////////////////////
// CASE 5 - create interpolated function from given data
void Readme_defining_functions_case_5_usage()
{
    // vectors containing values for interpolation
    Vector<double> x_val(100), y_val(100);

    LinearInterpRealFunc    f_linear(x_val, y_val);
    PolynomInterpRealFunc   f_polynom(x_val, y_val, 3);
    RationalInterpRealFunc  f_rational(x_val, y_val, 3);
    SplineInterpRealFunc    f_spline(x_val, y_val);
    BaryRatInterpRealFunc   f_baryrat(x_val, y_val, 3);
}

void Readme_defining_functions()
{
    Readme_defining_functions_case_1_usage();
    Readme_defining_functions_case_2_usage();
    Readme_defining_functions_case_3_usage();
    Readme_defining_functions_case_4_usage();
    Readme_defining_functions_case_5_usage();
}
