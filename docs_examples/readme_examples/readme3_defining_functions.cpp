#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector.h"
#include "base/Matrix.h"

#include "core/Function.h"
#include "core/InterpolatedFunction.h"
#include "core/Serializer.h"
#endif

#include "../test_data/real_functions_test_bed.h"
#include "../test_data/scalar_functions_test_bed.h"
#include "../test_data/vector_functions_test_bed.h"
#include "../test_data/parametric_curves_test_bed.h"

using namespace MML;

/////////////////////////////////////////////////////////////////////////////////////
// CASE 1 - standalone function providing calculation of a function
Real Readme_functions_TestFunc(Real x) { 
    return sin(x)*(1 + x*x / 2); 
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
    RealFunction f2{[](Real x) { return (Real) sin(x)*(1 + x*x / 2); } };

    // creating directly different types of functions
    ScalarFunction<3>       fScalar([](const VectorN<Real, 3> &x) { return x[0]; });
    VectorFunction<3>       fVector([](const VectorN<Real, 3> &x) { return VectorN<Real, 3>{0, x[0] * x[1], 0}; });
    VectorFunctionNM<2, 3>  fVectorNM([](const VectorN<Real, 2> &x) { return VectorN<Real, 3>{0, x[0] * x[1], 0}; });
    ParametricCurve<3>      paramCurve([](Real x) { return VectorN<Real, 3>{x, 2 * x, 3 * x}; });
    ParametricSurface<3>    paramSurface([](Real x, Real y) { return VectorN<Real, 3>{x * y, 2 * x * y, 3 * x}; });  

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
        Real operator()(Real x ) const { 
            return 1;     /* calculation using member variables */ 
        }
};
// Option 2 - inherit class from IRealFunction interface, and use the object itself as RealFunction
class ClassProvidingFuncToDerive2 : public IRealFunction {
    public:
        Real operator()(Real x ) const { 
            return 1.0;     /* calculation using member variables */ 
        }
};
void Readme_defining_functions_case_3_usage()
{
    ClassProvidingFuncToDerive   obj1;
    RealFunctionFromStdFunc f1(std::function<Real(Real)>{obj1});
    
    ClassProvidingFuncToDerive2   f2;     // usable RealFunction object (can be derived, integrated, ...)
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

    Real operator()(Real x ) const {
        return 1;     /* calculation using _ref */ 
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
    const int NumInterpPnt = 12;
    const Real x1 = 0, x2 = 10.0;

    // we will use this as test func
    RealFunction test_func{[](Real x) { return (Real) sin(x)*(1 + x*x / 2); } };
    
    // and using our helper, available for all real functions, create data for interpolation
    Vector<Real> x_val(NumInterpPnt), y_val(NumInterpPnt);
    test_func.GetValues(x1, x2, NumInterpPnt, x_val, y_val);

    // these are the ways we can interpolate Real function
    LinearInterpRealFunc    f_linear(x_val, y_val);
    PolynomInterpRealFunc   f_polynom(x_val, y_val, 3);
    SplineInterpRealFunc    f_spline(x_val, y_val);
    BaryRatInterpRealFunc   f_baryrat(x_val, y_val, 3);

    // situation - we need different number of points for different functions
    Serializer::SaveRealFuncEquallySpacedDetailed(test_func, "readme_interp_test_func", x1, x2, 100, "..\\..\\results\\readme_interp_test_func.txt");
    Serializer::SaveRealFuncEquallySpacedDetailed(f_linear, "readme_interp_linear_5_pnt", x1, x2, 500, "..\\..\\results\\readme_interp_linear_5_pnt.txt");
    Serializer::SaveRealFuncEquallySpacedDetailed(f_polynom, "readme_interp_polynom_5_pnt", x1, x2, 100, "..\\..\\results\\readme_interp_polynom_5_pnt.txt");
    Serializer::SaveRealFuncEquallySpacedDetailed, (x1, x2, 100, "..\\..\\results\\readme_interp_spline_5_pnt.txt");
    Serializer::SaveRealFuncEquallySpacedDetailed(f_baryrat, "readme_interp_baryrat_5_pnt", x1, x2, 100, "..\\..\\results\\readme_interp_baryrat_5_pnt.txt");

    const char *cmd = "..\\..\\tools\\visualizers\\real_function_visualizer\\MML_RealFunctionVisualizer.exe"
                        " ..\\..\\results\\readme_interp_test_func.txt"
                        " ..\\..\\results\\readme_interp_linear_5_pnt.txt"
                        " ..\\..\\results\\readme_interp_polynom_5_pnt.txt"
                        " ..\\..\\results\\readme_interp_spline_5_pnt.txt"
                        " ..\\..\\results\\readme_interp_baryrat_5_pnt.txt";
    std::system(cmd);

    // TODO 0.9 - HIGH parametric curve interpolation
}

void Readme_defining_functions()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                  README - defining functions                  ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Readme_defining_functions_case_1_usage();
    Readme_defining_functions_case_2_usage();
    Readme_defining_functions_case_3_usage();
    Readme_defining_functions_case_4_usage();
    Readme_defining_functions_case_5_usage();
}
