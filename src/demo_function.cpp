#include <iostream>
#include <iomanip>
#include <cmath>

#include "basic_types/Function.h"

double TestFunc1(double x)
{
    return x;
}

void Demo_Function()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                          FUNCTION                             ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    std::function<double(double)> f(TestFunc1);
    
    MML::RealFunctionFromStdFunc func(f);    
}