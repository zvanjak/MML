// TODO - BIG, EMPTY!!! - dovršiti demo
#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "basic_types/DiracDeltaFunction.h"
#endif

using namespace MML;


void Demo_Dirac_function()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                       DIRAC FUNCTION                          ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    DiracExp dirac1(10);

    double v = 3.0 * dirac1(0.0);
}