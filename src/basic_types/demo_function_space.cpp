
#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/FunctionSpace.h"
#endif

using namespace MML;

// TODO - BIG, finish Demo function space
void Demo_Function_Space()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                       FUNCTION SPACE                          ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    HermitianFunctionSpace5 hermite5(0.0, 1.0);
    LegendreFunctionSpace5  legendre5;

    RealFunction f([](Real x){ return x * sin(x);});

    auto rep = hermite5.getRepr(f);
    std::cout << "Representation of function x * sin(x)\n";
    rep.Print(std::cout, 10, 6);

    Real x = 0.3;
    auto v = rep[0] * hermite5.getBasisFunc(0)(x) + rep[1] * hermite5.getBasisFunc(1)(x) + rep[2] * hermite5.getBasisFunc(2)(x) + rep[3] * hermite5.getBasisFunc(3)(x) + rep[4] * hermite5.getBasisFunc(4)(x);
    std::cout << "\nReconstructed function from its representation\n";
    std::cout << "f(" << x << ") = " << x * sin(x)  << std::endl;
    std::cout << "f(" << x << ") = " << v << std::endl;

}