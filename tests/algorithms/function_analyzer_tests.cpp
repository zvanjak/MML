#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"

#include "core/Function.h"
#include "core/InterpolatedFunction.h"

#include "algorithms/FunctionAnalyzer.h"
#endif

using namespace MML;

double test_func_123(double x)
{
    return 9 - (x-3)*(x-3);
}

// TODO - finalize these tests
// POUZDANA METRIKA RAZLIKE DVIJE FUNKCIJE!!!!!

TEST_CASE("Test_Function_Analyzer_FuncDiff", "[simple]") 
{
    Vector<Real> x{0.0, 3.0, 6.0};
    Vector<Real> y{test_func_123(0.0), test_func_123(3.0), test_func_123(6.0)};

    LinearInterpRealFunc myfunc(x,y);
    RealFunction test(test_func_123);
    
    double triangleArea = 6 * test_func_123(3.0) / 2;
    double parabolaArea = Integration::IntegrateTrap(test, 0.0, 6.0);

    REQUIRE( FunctionComparer::getIntegratedDiff(myfunc, test, 0.0, 6.0) == Approx(-9.0).epsilon(1e-6) );
    REQUIRE( FunctionComparer::getIntegratedDiff(myfunc, test, 0.0, 6.0) != Approx(-9.0).epsilon(1e-7) );

}
