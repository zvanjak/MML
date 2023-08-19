#include <iostream>
#include <iomanip>
#include <cmath>

#include <vector>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "basic_types/InterpolatedFunction.h"
#endif

using namespace MML;

double TestFunc2(double x)
{
    return sin(x);
}

void Demo_Interpolated_Function()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                     INTERPOLATED FUNCTION                     ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // std::vector<double> values(101);

    // for(int i=0; i<=100; i++ )
    //     values[i] = sin(i / 10.0);

    // MML::TabulatedValues1DEqualSpacing tabValues{0.0, 10.0, values};


}