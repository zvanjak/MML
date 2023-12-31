#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Vector.h"
#include "core/Matrix.h"

#include "core/Function.h"

#endif

using namespace MML;

double Readme_functions_TestFunc(double x) 
{ 
    return sin(x)*(1.0 + 0.5*x*x); 
}
void Readme_functions()
{
    // creating a function object from a already existing (standalone) function
    RealFunction f1(Readme_functions_TestFunc);

    // or creating a function object directly
    RealFunction f2{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };
}
