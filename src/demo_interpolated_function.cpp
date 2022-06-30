#include <iostream>
#include <iomanip>
#include <cmath>

#include <vector>

#include "basic_types/InterpolatedFunction.h"

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

    std::vector<double> values(101);

    for(int i=0; i<=100; i++ )
        values[i] = sin(i / 10.0);

    MML::TabulatedValues1DEqualSpacing tabValues{0.0, 10.0, values};


}